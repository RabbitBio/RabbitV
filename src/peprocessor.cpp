#include "peprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "adaptertrimmer.h"
#include "basecorrector.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "polyx.h"

// modified
#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include <string>
#include "io/DataQueue.h"
#include "io/Formater.h"
typedef mash::core::TDataQueue<mash::fq::FastqDataPairChunk>FqPairChunkQueue;
typedef mash::core::TDataQueue<mash::fq::FastqDataChunk>FqChunkQueue;
// modified over

PairEndProcessor::PairEndProcessor(Options* opt){
    readNum = 0;
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream1 = NULL;
    mZipFile1 = NULL;
    mOutStream2 = NULL;
    mZipFile2 = NULL;
    mUmiProcessor = new UmiProcessor(opt);

    int isizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new long[isizeBufLen];
    memset(mInsertSizeHist, 0, sizeof(long)*isizeBufLen);
    mLeftWriter =  NULL;
    mRightWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }

    mVirusDetector = new VirusDetector(opt);
}

PairEndProcessor::~PairEndProcessor() {
    delete mInsertSizeHist;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    if(mVirusDetector) {
        delete mVirusDetector;
        mVirusDetector = NULL;
    }
}

void PairEndProcessor::initOutput() {
    if(mOptions->out1.empty())
        return;
    
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    if(!mOptions->out2.empty())
        mRightWriter = new WriterThread(mOptions, mOptions->out2);
}

void PairEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mRightWriter) {
        delete mRightWriter;
        mRightWriter = NULL;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
}


bool PairEndProcessor::process(){
    initOutput();

    initPackRepository();
    // modified
    mash::fq::FastqDataPool* fastqPool = new mash::fq::FastqDataPool(256, 1 << 22);
    FqPairChunkQueue queue1(128, 1);
    FqChunkQueue queue2(128,1);
    // modified over

    std::thread producer(std::bind(&PairEndProcessor::producerTask, this,fastqPool,std::ref(queue1),std::ref(queue2)));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, t, true);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this, configs[t],fastqPool,std::ref(queue1),std::ref(queue2)));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mLeftWriter));
    if(mRightWriter)
        rightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mRightWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(leftWriterThread)
        leftWriterThread->join();
    if(rightWriterThread)
        rightWriterThread->join();

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and filter results
    vector<Stats*> preStats1;
    vector<Stats*> postStats1;
    vector<Stats*> preStats2;
    vector<Stats*> postStats2;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats1.push_back(configs[t]->getPreStats1());
        postStats1.push_back(configs[t]->getPostStats1());
        preStats2.push_back(configs[t]->getPreStats2());
        postStats2.push_back(configs[t]->getPostStats2());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats1 = Stats::merge(preStats1);
    Stats* finalPostStats1 = Stats::merge(postStats1);
    Stats* finalPreStats2 = Stats::merge(preStats2);
    Stats* finalPostStats2 = Stats::merge(postStats2);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

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
        cerr << "Duplication rate: " << dupRate * 100.0 << "%" << endl;
    }

    // insert size distribution
    int peakInsertSize = getPeakInsertSize();
    cerr << endl;
    cerr << "Insert size peak (evaluated by paired-end reads): " << peakInsertSize << endl;

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(mVirusDetector, finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.setInsertHist(mInsertSizeHist, peakInsertSize);
    hr.report(mVirusDetector, finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats1;
    delete finalPostStats1;
    delete finalPreStats2;
    delete finalPostStats2;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;
    if(rightWriterThread)
        delete rightWriterThread;

    closeOutput();

    return true;
}

int PairEndProcessor::getPeakInsertSize() {
    int peak = 0;
    long maxCount = -1;
    for(int i=0; i<mOptions->insertSizeMax; i++) {
        if(mInsertSizeHist[i] > maxCount) {
            peak = i;
            maxCount = mInsertSizeHist[i];
        }
    }
    return peak;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config){
    string outstr1;
    string outstr2;
    string singleOutput;
    int readPassed = 0;
    int mergedCount = 0;
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->mLeft;
        Read* or2 = pair->mRight;

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statPair(or1, or2);

        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1, or2);

        // trim in head and tail, and apply quality cut in sliding window
        int frontTrimmed1 = 0;
        int frontTrimmed2 = 0;
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed1);
        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2, frontTrimmed2);

        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, r2, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }
        bool isizeEvaluated = false;
        if(r1 != NULL && r2!=NULL && (mOptions->adapter.enabled || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            // we only use thread 0 to evaluae ISIZE
            if(config->getThreadId() == 0) {
                statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
                isizeEvaluated = true;
            }
            if(mOptions->correction.enabled) {
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            if(mOptions->adapter.enabled) {
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov, frontTrimmed1, frontTrimmed2);
                bool trimmed1 = trimmed;
                bool trimmed2 = trimmed;
                if(!trimmed){
                    if(mOptions->adapter.hasSeqR1)
                        trimmed1 = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
                    if(mOptions->adapter.hasSeqR2)
                        trimmed2 = AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.sequenceR2, true);
                }
                if(mOptions->adapter.hasFasta) {
                    AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, !trimmed1);
                    AdapterTrimmer::trimByMultiSequences(r2, config->getFilterResult(), mOptions->adapter.seqsInFasta, true, !trimmed2);
                }
            }
        }

        if(config->getThreadId() == 0 && !isizeEvaluated && r1 != NULL && r2!=NULL) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit/100.0);
            statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
            isizeEvaluated = true;
        }

        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, r2, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL && r2!=NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
            if( mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length())
                r2->resize(mOptions->trim.maxLen2);
        }


        int result1 = mFilter->passFilter(r1);
        int result2 = mFilter->passFilter(r2);

        config->addFilterResult(max(result1, result2), 2);

        if( r1 != NULL &&  result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER ) {

            bool found = false;
            found |= mVirusDetector->detect(r1);
            found |= mVirusDetector->detect(r2);
            
            if(found) {
                if(mOptions->outputToSTDOUT) {
                    singleOutput += r1->toString() + r2->toString();
                } else {
                    outstr1 += r1->toString();
                    outstr2 += r2->toString();
                }
            }

            // stats the read after filtering
            config->getPostStats1()->statRead(r1);
            config->getPostStats2()->statRead(r2);

            readPassed++;
        }

        delete pair;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            delete r1;
        // if no trimming applied, r1 should be identical to or1
        if(r2 != or2 && r2 != NULL)
            delete r2;
    }

    // normal output by left/right writer thread
    if(mRightWriter && mLeftWriter && (!outstr1.empty() || !outstr2.empty())) {
        // write PE
        char* ldata = new char[outstr1.size()];
        memcpy(ldata, outstr1.c_str(), outstr1.size());
        mLeftWriter->input(ldata, outstr1.size());

        char* rdata = new char[outstr2.size()];
        memcpy(rdata, outstr2.c_str(), outstr2.size());
        mRightWriter->input(rdata, outstr2.size());
    } else if(mLeftWriter && !singleOutput.empty()) {
        // write singleOutput
        char* ldata = new char[singleOutput.size()];
        memcpy(ldata, singleOutput.c_str(), singleOutput.size());
        mLeftWriter->input(ldata, singleOutput.size());
    }

    mOutputMtx.unlock();

    config->markProcessed(pack->count);

    delete pack->data;
    delete pack;

    return true;
}
    
void PairEndProcessor::statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1, int frontTrimmed2) {
    int isize = mOptions->insertSizeMax;
    if(ov.overlapped) {
        if(ov.offset > 0)
            isize = r1->length() + r2->length() - ov.overlap_len + frontTrimmed1 + frontTrimmed2;
        else
            isize = ov.overlap_len + frontTrimmed1 + frontTrimmed2;
    }

    if(isize > mOptions->insertSizeMax)
        isize = mOptions->insertSizeMax;

    mInsertSizeHist[isize]++;
}

bool PairEndProcessor::processRead(Read* r, ReadPair* originalPair, bool reversed) {
    // do something here
    return true;
}

void PairEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPairPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPairPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    
}

void PairEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void PairEndProcessor::producePack(ReadPairPack* pack){
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

// consumePack modified
/*
void PairEndProcessor::consumePack(ThreadConfig* config){
    ReadPairPack* data;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    } * /

    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if(mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;* /
    mInputMtx.unlock();
    //mRepo.readPos++;

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();

    processPairEnd(data, config);

}
*/
void PairEndProcessor::consumePack(ThreadConfig* config, ReadPairPack* data) {
    processPairEnd(data, config);
}
// consumePack modified over
// modified
/*
void PairEndProcessor::producerTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    ReadPair** data = new ReadPair*[PACK_SIZE];
    memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput);
    int count=0;
    bool needToBreak = false;
    while(true){
        ReadPair* read = reader.read();
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(!read || needToBreak){
            // the last pack
            ReadPairPack* pack = new ReadPairPack;
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
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M read pairs";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new ReadPair*[PACK_SIZE];
            memset(data, 0, sizeof(ReadPair*)*PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                slept++;
                usleep(1000);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while( (mLeftWriter && mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) || (mRightWriter && mRightWriter->bufferLength() > PACK_IN_MEM_LIMIT) ){
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
                    reader.mLeft->getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }* /
        }
    }

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete[] data;
}
*/
void PairEndProcessor::producerTask(mash::fq::FastqDataPool* fastqPool, FqPairChunkQueue& dq,FqChunkQueue& dq2)
{
    if (mOptions->verbose)
        loginfo("start to load data");
    mash::fq::FastqFileReader* fqFileReader;
    bool isZipped = false;
    if (ends_with(mOptions->in1, ".gz")) {
        isZipped = true;
    }
    fqFileReader = new mash::fq::FastqFileReader(mOptions->in1, *fastqPool, mOptions->in2, isZipped);
    int n_chunks = 0;
    if(mOptions->interleavedInput){
        while(true){
            mash::fq::FastqChunk * fqChunk2 = new mash::fq::FastqChunk;
            fqChunk2->chunk = fqFileReader->readNextPairChunkInterleaved();
            if(fqChunk2->chunk == NULL)
                break;
            n_chunks++;
            dq2.Push(n_chunks,fqChunk2->chunk);
        }
        dq2.SetCompleted();
    }else{
        while (true) {
            mash::fq::FastqPairChunk* fqChunk = new mash::fq::FastqPairChunk;
            fqChunk->chunk = fqFileReader->readNextPairChunk();
            if (fqChunk->chunk == NULL)
                break;
            n_chunks++;
            dq.Push(n_chunks, fqChunk->chunk);
        }
        dq.SetCompleted();
     }

    delete fqFileReader;
    std::cout << "file " << mOptions->in1 << " has " << n_chunks << " chunks" << std::endl;
    mProduceFinished = true;
    if (mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
}
// modified over

// modified
/*
void PairEndProcessor::consumerTask(ThreadConfig* config)
{
    while(true) {
        if(config->canBeStopped()){
            mFinishedThreads++;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos) {
            if(mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            mFinishedThreads++;
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        }
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mRightWriter)
            mRightWriter->setInputCompleted();
    }
    
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}
*/
void PairEndProcessor::consumerTask(ThreadConfig* config, mash::fq::FastqDataPool* fastqPool, FqPairChunkQueue & dq,FqChunkQueue& dq2)
{
    // TODO when file is interleaved
    bool needToBreak = false;
    mash::int64 id = 0;
    std::vector<neoReference> left_data;
    std::vector<neoReference> right_data;
    left_data.reserve(10000);
    right_data.reserve(10000);
    mash::fq::FastqPairChunk* fqChunk = new mash::fq::FastqPairChunk;
    mash::fq::FastqChunk* left_fqChunk = new mash::fq::FastqChunk;
    mash::fq::FastqChunk* right_fqChunk = new mash::fq::FastqChunk;
    mash::int64 left_start = 0, right_start = 0;
    mash::int64 left_count = 0, right_count = 0;
    mash::int64 left_end  = left_start + left_count;
    ReadPair** data = new ReadPair * [PACK_SIZE];
    memset(data, 0, sizeof(ReadPair*) * PACK_SIZE);
    int count = 0; // 当前data中readPair的数量
    // interleaved
    mash::fq::FastqChunk* fqChunk2 = new mash::fq::FastqChunk;
    if(mOptions->interleavedInput){
        while(dq2.Pop(id,fqChunk2->chunk)){
            if(config->canBeStopped()){
                break;
            }
            left_count = mash::fq::chunkFormat(fqChunk2,left_data,true);
            left_end = left_start + left_count;
            std::string name1 = std::string((char*)left_data[left_start].base + left_data[left_start].pname, left_data[left_start].lname);
            std::string name2 = std::string((char*)left_data[left_start+1].base + left_data[left_start+1].pname, left_data[left_start+1].lname);
            // 靠靠read靠靠靠pair
            //if(name1.substr(0,name1.find(' ')) != name2.substr(0,name2.find(' '))){
            if(std::stol(name1.substr(16,name1.find(' '))) % 2 == 0){
                //printf("name1: %s\n",(name1.substr(0,name1.find(' '))));
                //std::cout<<"name :"<<name1.substr(0,16)<<std::endl;
                left_start++;
            }

            for(mash::int64 start = left_start; start < left_end - 1; start++){
                std::string name1 = std::string((char*)left_data[start].base + left_data[start].pname, left_data[start].lname);
                std::string seq1 = std::string((char*)left_data[start].base + left_data[start].pseq, left_data[start].lseq);
                std::string strand1 = std::string((char*)left_data[start].base + left_data[start].pstrand, left_data[start].lstrand);
                std::string quality1 = std::string((char*)left_data[start].base + left_data[start].pqual, left_data[start].lqual);
                Read* l = new Read(name1, seq1, strand1, quality1, mOptions->phred64);
                start++;
                std::string name2 = std::string((char*)left_data[start+1].base + left_data[start+1].pname, left_data[start+1].lname);
                std::string seq2 = std::string((char*)left_data[start].base + left_data[start].pseq, left_data[start].lseq);
                std::string strand2 = std::string((char*)left_data[start].base + left_data[start].pstrand, left_data[start].lstrand);
                std::string quality2 = std::string((char*)left_data[start].base + left_data[start].pqual, left_data[start].lqual);
                Read* r = new Read(name2, seq2, strand2, quality2, mOptions->phred64);
                ReadPair* read = new ReadPair(l, r);
                data[count] = read;
                count++;
                // configured to process only first N reads
                if (mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess) {
                    needToBreak = true;
                }
                if (count == PACK_SIZE || needToBreak) {
                    ReadPairPack* pack = new ReadPairPack;
                    pack->data = data;
                    pack->count = count;
                    readNum += count;
                    // consumPack
                    consumePack(config,pack);
                    data = new ReadPair * [PACK_SIZE];
                    memset(data, 0, sizeof(ReadPair*) * PACK_SIZE);
                    count = 0;
                }

            }// for
            left_start  = left_end;
            fastqPool->Release(fqChunk2->chunk);
        } // while
    }else{
        // not interleaved
        while (dq.Pop(id, fqChunk->chunk)) {
            if (config->canBeStopped()) {
                // mFinishedThreads++;
                break;
            }
            left_fqChunk->chunk = fqChunk->chunk->left_part;
            right_fqChunk->chunk = fqChunk->chunk->right_part;
            left_count = mash::fq::chunkFormat(left_fqChunk, left_data, true);
            right_count = mash::fq::chunkFormat(right_fqChunk, right_data, true);
            // 左右datachunk里面的read数量相等
            for (mash::int64 start = left_start; start < left_start + left_count; start++) {
                std::string name1 = std::string((char*)left_data[start].base + left_data[start].pname, left_data[start].lname);
                std::string seq1 = std::string((char*)left_data[start].base + left_data[start].pseq, left_data[start].lseq);
                std::string strand1 = std::string((char*)left_data[start].base + left_data[start].pstrand, left_data[start].lstrand);
                std::string quality1 = std::string((char*)left_data[start].base + left_data[start].pqual, left_data[start].lqual);
                Read* l = new Read(name1, seq1, strand1, quality1, mOptions->phred64);
                std::string name2 = std::string((char*)right_data[start].base + right_data[start].pname, right_data[start].lname);
                std::string seq2 = std::string((char*)right_data[start].base + right_data[start].pseq, right_data[start].lseq);
                std::string strand2 = std::string((char*)right_data[start].base + right_data[start].pstrand, right_data[start].lstrand);
                std::string quality2 = std::string((char*)right_data[start].base + right_data[start].pqual, right_data[start].lqual);
                Read* r = new Read(name2, seq2, strand2, quality2, mOptions->phred64);
                ReadPair* read = new ReadPair(l, r);
                data[count] = read;
                count++;
                // configured to process only first N reads
                if (mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess) {
                    needToBreak = true;
                }
                if (count == PACK_SIZE || needToBreak) {
                    ReadPairPack* pack = new ReadPairPack;
                    pack->data = data;
                    pack->count = count;
                    readNum += count;
                    // consumPack
                    consumePack(config,pack);
                    data = new ReadPair * [PACK_SIZE];
                    memset(data, 0, sizeof(ReadPair*) * PACK_SIZE);
                    count = 0;
                }

            } // for
    
            left_start += left_count;
            fastqPool->Release(fqChunk->chunk->left_part);
            fastqPool->Release(fqChunk->chunk->right_part);
        } // while
    }// interleaved else
    delete fqChunk2;
    delete left_fqChunk;
    delete right_fqChunk;
    // 从dataqueue中获取不到新的chunk的时候 判断之前data是否还存在没处理的数据
    if (count > 0) {
        ReadPairPack* pack = new ReadPairPack;
        pack->data = data;
        pack->count = count;
        consumePack(config,pack);
        readNum += count;
    }
    mFinishedThreads++;
    if (mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
        loginfo(msg);
    }
    if (mFinishedThreads == mOptions->thread) {
        if (mLeftWriter)
            mLeftWriter->setInputCompleted();
        if (mRightWriter)
            mRightWriter->setInputCompleted();
    }

    if (mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }

}

// modified over

void PairEndProcessor::writeTask(WriterThread* config)
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
