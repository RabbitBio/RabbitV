#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"
// modified
#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include <string>
#include "io/DataQueue.h"
#include "io/Formater.h"
typedef mash::core::TDataQueue<mash::fq::FastqDataChunk> FqChunkQueue;
// modified over

SingleEndProcessor::SingleEndProcessor(Options* opt){
    readNum = 0;
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream = NULL;
    mZipFile = NULL;
    mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter =  NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }

    mVirusDetector = new VirusDetector(opt);
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    if(mVirusDetector) {
        delete mVirusDetector;
        mVirusDetector = NULL;
    }
}

void SingleEndProcessor::initOutput() {
    if(mOptions->out1.empty())
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
}

bool SingleEndProcessor::process(){
    initOutput();

    initPackRepository();

    //modified 
    mash::fq::FastqDataPool* fastqPool = new mash::fq::FastqDataPool(256, 1 << 22);
    FqChunkQueue queue1(128, 1);

    //modified over
    std::thread producer(std::bind(&SingleEndProcessor::producerTask, this,fastqPool,std::ref(queue1)));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, t, false);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t],fastqPool,std::ref(queue1)));
    }

    std::thread* leftWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(leftWriterThread)
        leftWriterThread->join();

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    // read filter results to the first thread's
    for(int t=1; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
    }

    mVirusDetector->report();

    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    }

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(mVirusDetector, finalFilterResult, finalPreStats, finalPostStats);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(mVirusDetector, finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;

    closeOutput();

    return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){
    string outstr;
    string failedOut;
    int readPassed = 0;
    for(int p=0;p<pack->count;p++){

        // original read1
        Read* or1 = pack->data[p];

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statRead(or1);

        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1);

        int frontTrimmed = 0;
        // trim in head and tail, and apply quality cut in sliding window
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed);

        if(r1 != NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }

        if(r1 != NULL && mOptions->adapter.enabled){
            bool trimmed = false;
            if(mOptions->adapter.hasSeqR1)
                trimmed = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
            bool incTrimmedCounter = !trimmed;
            if(mOptions->adapter.hasFasta) {
                AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, incTrimmedCounter);
            }
        }

        if(r1 != NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
        }

        int result = mFilter->passFilter(r1);

        config->addFilterResult(result, 1);

        if( r1 != NULL &&  result == PASS_FILTER) {

            bool found = mVirusDetector->detect(r1);

            if(found)
                outstr += r1->toString();

            // stats the read after filtering
            config->getPostStats1()->statRead(r1);
            readPassed++;
        }

        delete or1;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            delete r1;
    }
    // if splitting output, then no lock is need since different threads write different files
    if(mLeftWriter) {
        char* ldata = new char[outstr.size()];
        memcpy(ldata, outstr.c_str(), outstr.size());
        mLeftWriter->input(ldata, outstr.size());
    }

    mOutputMtx.lock();
    if(mOptions->outputToSTDOUT) {
        fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    }
/*
    if(mLeftWriter) {
        char* ldata = new char[outstr.size()];
        memcpy(ldata, outstr.c_str(), outstr.size());
        mLeftWriter->input(ldata, outstr.size());
    }
*/
    mOutputMtx.unlock();

    config->markProcessed(pack->count);

    delete pack->data;
    delete pack;

    return true;
}

void SingleEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    //mRepo.readCounter = 0;
    
}

void SingleEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack){
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    /*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        //mRepo.repoNotFull.wait(lock);
    }*/

    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;

    /*if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;*/

    //mRepo.repoNotEmpty.notify_all();
    //lock.unlock();
}

void SingleEndProcessor::consumePack(ThreadConfig* config,ReadPack* pack){
    ReadPack* data;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    }*/

    // mmodified
    //mInputMtx.lock();
    //while(mRepo.writePos <= mRepo.readPos) {
    //    usleep(1000);
    //    if(mProduceFinished) {
    //        mInputMtx.unlock();
    //        return;
    //    }
    //}
    //data = mRepo.packBuffer[mRepo.readPos]; //读一个pack
    //mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    //mInputMtx.unlock();

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();
    data = pack;

    processSingleEnd(data, config);

}

void SingleEndProcessor::producerTask(mash::fq::FastqDataPool* fastqPool, FqChunkQueue& dq)
{
    if(mOptions->verbose)
        loginfo("start to load data"); //debug information
    //long lastReported = 0;
    //int slept = 0; 
    //long readNum = 0; // 之前累积读取的read的个数
    //bool splitSizeReEvaluated = false;
    //Read** data = new Read*[PACK_SIZE];
    //memset(data, 0, sizeof(Read*)*PACK_SIZE);
    //FastqReader reader(mOptions->in1, true, mOptions->phred64); // 
    
    // modified
    mash::fq::FastqFileReader* fqFileReader;
    // 是否是压缩文件
    bool isZipped = false;
    if (ends_with(mOptions->in1, ".gz")) {
        isZipped = true;
    }
    fqFileReader = new mash::fq::FastqFileReader(mOptions->in1, (*fastqPool), "",isZipped);
    int n_chunks = 0;
    int line_sum = 0;
    while (true) {
        mash::fq::FastqChunk* fqChunk = new mash::fq::FastqChunk;
        fqChunk->chunk = fqFileReader->readNextChunk();
        if (fqChunk->chunk == NULL)
            break;
        n_chunks++;
       // std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqChunk->chunk);
    }
    dq.SetCompleted();
    delete fqFileReader;
    std::cout << "file " << mOptions->in1 << " has " << n_chunks << " chunks" << std::endl;
    // modified over
    
    /*
    int count=0; // data or pack中read的个数
    bool needToBreak = false;
    while(true){
        Read* read = reader.read(); // read()方法返回一个read对象指针
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(!read || needToBreak){
            // the last pack
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads 只读入前readToProcess条read
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        // debug时每百万条read后打印一次信息
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M reads";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            // 读的过快时 等待consumer 为什么要等？？？
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                //cerr<<"sleep"<<endl;
                slept++;
                usleep(100);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            // 500 packs
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while(mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) {
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
            // re-evaluate split size
            // TODO: following codes are commented since it may cause threading related conflicts in some systems
            /*if(mOptions->split.needEvaluation && !splitSizeReEvaluated && readNum >= mOptions->split.size) {
                splitSizeReEvaluated = true;
                // greater than the initial evaluation
                if(readNum >= 1024*1024) {
                    size_t bytesRead;
                    size_t bytesTotal;
                    reader.getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }* /
        }
    }//读取read完成
    */
    
    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();

    // if the last data initialized is not used, free it
    //if(data != NULL)
    //    delete[] data;
}

void SingleEndProcessor::consumerTask(ThreadConfig* config,mash::fq::FastqDataPool* fastqPool, FqChunkQueue& dq)
{
    //while(true) {
    //    if(config->canBeStopped()){
    //        mFinishedThreads++;
    //        break;
    //    }

        // modified
        mash::int64 seq_count_start = 0;
        mash::int64 seq_count = 0;
        mash::int64 id = 0;
        std::vector<neoReference> data;
        mash::fq::FastqChunk* fqChunk = new mash::fq::FastqChunk;
        data.reserve(10000);
        Read** data2 = new Read*[PACK_SIZE];
        memset(data2, 0, sizeof(Read*)*PACK_SIZE);
        int count = 0; // data or pack中read的个数
        bool needToBreak = false;
        while (dq.Pop(id, fqChunk->chunk)) {
            seq_count = mash::fq::chunkFormat(fqChunk, data, true);
            //Read** data2 = new Read*[seq_count];
            //memset(data2, 0, sizeof(Read*)*seq_count);
            // 将该FastqDataChunk进行format 转化成pack的形式 并存储到mRepo中
            for (mash::int64 start = seq_count_start; start < seq_count_start + seq_count; start++) {
                std::string name = std::string((char*)data[start].base + data[start].pname, data[start].lname);
                std::string seq = std::string((char*)data[start].base + data[start].pseq, data[start].lseq);
                std::string strand = std::string((char*)data[start].base + data[start].pstrand, data[start].lstrand);
                std::string quality = std::string((char*)data[start].base + data[start].pqual, data[start].lqual);
                Read* read = new Read(name, seq, strand, quality, mOptions->phred64);
                data2[count] = read;
                count++;
                // configured to process only first N reads 只读入前readToProcess条read
                // TODO
                if (mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess) {
                    needToBreak = true;
                }
                // a full pack
                if (count == PACK_SIZE || needToBreak) {
                    ReadPack* pack = new ReadPack;
                    pack->data = data2;
                    pack->count = count;
                    //producePack(pack);
                    consumePack(config,pack);
                    //re-initialize data for next pack
                    data2 = new Read * [PACK_SIZE];
                    memset(data2, 0, sizeof(Read*) * PACK_SIZE);
                    readNum += count;
                    // reset count to 0
                    count = 0;
                }
            }
            seq_count_start += seq_count;
            fastqPool->Release(fqChunk->chunk);
        }
        if (count > 0) {
            ReadPack* pack = new ReadPack;
            pack->data = data2;
            pack->count = count;
            consumePack(config, pack);
            readNum += count; 
        }
        //// debug时每百万条read后打印一次信息
        //if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
        //    lastReported = count + readNum;
        //    string msg = "loaded " + to_string((lastReported/1000000)) + "M reads";
        //    loginfo(msg);
        //}
        
        //modified over

        //while(mRepo.writePos <= mRepo.readPos) {
        //    if(mProduceFinished)
        //        break;
        //    usleep(1000); // 等待生产者产生read
        //}
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        //if(mProduceFinished && mRepo.writePos == mRepo.readPos){
        //    mFinishedThreads++;
        //    if(mOptions->verbose) {
        //        string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
        //        loginfo(msg);
        //    }
        //    //lock.unlock();
        //    break;
        //}
        //if(mProduceFinished){
        //    if(mOptions->verbose) {
        //        string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
        //        loginfo(msg);
        //    }
        //    consumePack(config);
        //    //lock.unlock();
        //} else {
        //    //lock.unlock();
        //    consumePack(config);
        //}
    //}
    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mFinishedThreads++;
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
        loginfo(msg);
    }
    //lock.unlock();
    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
    }

    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::writeTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
