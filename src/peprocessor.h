#ifndef PE_PROCESSOR_H
#define PE_PROCESSOR_H

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
#include "overlapanalysis.h"
#include "writerthread.h"
#include "duplicate.h"
#include "virusdetector.h"

// modified
#include "io/FastxStream.h"
#include "io/FastxChunk.h"
#include <string>
#include "io/DataQueue.h"
#include "io/Formater.h"
typedef mash::core::TDataQueue<mash::fq::FastqDataPairChunk> FqPairChunkQueue;
typedef mash::core::TDataQueue<mash::fq::FastqDataChunk> FqChunkQueue;
// modified over


using namespace std;

struct ReadPairPack {
    ReadPair** data;
    int count;
};

typedef struct ReadPairPack ReadPairPack;

struct ReadPairRepository {
    ReadPairPack** packBuffer;
    atomic_long readPos;
    atomic_long writePos;
    //std::mutex mtx;
    //std::mutex readCounterMtx;
    //std::condition_variable repoNotFull;
    //std::condition_variable repoNotEmpty;
};

typedef struct ReadPairRepository ReadPairRepository;

class PairEndProcessor{
public:
    PairEndProcessor(Options* opt);
    ~PairEndProcessor();
    bool process();

private:
    bool processPairEnd(ReadPairPack* pack, ThreadConfig* config);
    bool processRead(Read* r, ReadPair* originalRead, bool reversed);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPairPack* pack);
    // modified
    //void consumePack(ThreadConfig* config);
    void consumePack(ThreadConfig* config, ReadPairPack* data);
    //void producerTask();
    void producerTask(mash::fq::FastqDataPool* fastqPool, FqPairChunkQueue& dq,FqChunkQueue& dq2);
    // void consumerTask(ThreadConfig* config);
    void consumerTask(ThreadConfig* config,mash::fq::FastqDataPool* fastqPool, FqPairChunkQueue& dq,FqChunkQueue& dq2);
    // modified over
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1 = 0, int frontTrimmed2 = 0);
    int getPeakInsertSize();
    void writeTask(WriterThread* config);

private:
    ReadPairRepository mRepo;
    atomic_bool mProduceFinished;
    atomic_int mFinishedThreads;
    std::mutex mOutputMtx;
    std::mutex mInputMtx;
    Options* mOptions;
    Filter* mFilter;
    gzFile mZipFile1;
    gzFile mZipFile2;
    ofstream* mOutStream1;
    ofstream* mOutStream2;
    UmiProcessor* mUmiProcessor;
    long* mInsertSizeHist;
    WriterThread* mLeftWriter;
    WriterThread* mRightWriter;
    Duplicate* mDuplicate;
    VirusDetector* mVirusDetector;
    atomic_long readNum;
};


#endif
