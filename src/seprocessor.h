#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>
#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "duplicate.h"
#include "virusdetector.h"
#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include "io/DataQueue.h"
#include "io/Formater.h"

#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include "io/Formater.h"
#include "io/Reference.h"

using namespace std;

typedef rabbit::fq::FastqFileReader FastqFileReader;
typedef rabbit::fq::FastqDataQueue FastqDataQueue;
typedef rabbit::fq::FastqDataPool FastqDataPool;
typedef rabbit::fq::FastqDataChunk FastqDataChunk;

struct ReadPack {
    Read** data;
    int count;
};

typedef struct ReadPack ReadPack;
//typedef mash::core::TDataQueue<mash::fq::FastqDataChunk> FqChunkQueue;

struct ReadRepository {
    ReadPack** packBuffer;
    atomic_long readPos;
    atomic_long writePos;
    //std::mutex mtx;
    //std::mutex readCounterMtx;
    //std::condition_variable repoNotFull;
    //std::condition_variable repoNotEmpty;
};

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor{
public:
    SingleEndProcessor(Options* opt);
    ~SingleEndProcessor();
    bool process();

private:
    bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPack* pack);
    void consumePack(ThreadConfig* config, ReadPack* pack);
    //void producerTask(mash::fq::FastqDataPool* fastqPool, FqChunkQueue& dq);
    void producerTask(FastqDataPool&, FastqDataQueue&);
    //void consumerTask(ThreadConfig* config,mash::fq::FastqDataPool* fastqPool,FqChunkQueue& dq);
    //
    void consumerTask(ThreadConfig* config, FastqDataPool&, FastqDataQueue&);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void writeTask(WriterThread* config);

private:
    Options* mOptions;
    ReadRepository mRepo;
    atomic_bool mProduceFinished;
    atomic_int mFinishedThreads;
    std::mutex mInputMtx;
    std::mutex mOutputMtx;
    Filter* mFilter;
    gzFile mZipFile;
    ofstream* mOutStream;
    UmiProcessor* mUmiProcessor;
    WriterThread* mLeftWriter;
    Duplicate* mDuplicate;
    VirusDetector* mVirusDetector;
    atomic_long readNum;
};


#endif
