#ifndef ALLKMER_H
#define ALLKMER_H

// includes
#include "common.h"
#include <vector>
//#include <unordered_map>
#include <map>
#include "fastareader.h"
#include "options.h"
#include "zlib/zlib.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <mutex>
#include <deque>
#include <thread>
#include "./io/DataQueue.h"

//change
#include <atomic>
#include "omp.h"

class padding_atomic_lock
{
    public:
        std::atomic_flag at_lock;
        char padding_arr[63];
};
//change

#define  MTX_COUNT 100
#define COLLISION_FLAG 0xFFFFFFFF

using namespace std;

typedef rabbit::core::TDataQueue<std::tuple<uint32, uint64_t*, uint64> > data_queue_t;

class KCResult {
public:
    string mName;
    uint64  mHit;
    int mMedianHit;
    double mMeanHit;
    double mCoverage;
    int mKmerCount;
    int mUniqueReads;
};

struct KCHit{
  uint64 mKey64;
  uint32 mID;
  uint32 mHit;
};

struct KC_t{
  uint32 mID;
  uint32 mHit;
};

class KmerCollection
{
public:
    bool readBin = false;
public:
    KmerCollection(string filename, Options* opt);
    ~KmerCollection();
    void init();
    void report();
    void reportJSON(ofstream& ofs);
    void reportHTML(ofstream& ofs);
    uint32 add(uint64 kmer64);
    uint32 add_bin(uint64 kmer64);
    void addGenomeRead(uint32 genomeID);
    void mul_thread_init();
    void readAllBin();

    uint32 packIdCount(uint32 id, uint32 count);
    void unpackIdCount(uint32 data,uint32& id, uint32& count);
    void stat();

private:
    bool getLine(char* line, int maxLine);
    uint64 makeHash(uint64 key);
    bool eof();
    void makeBitAndMask();
    bool isHighConfidence(KCResult kcr);
private:
    Options* mOptions;
    vector<string> mNames;
    vector<uint64> mHits;
    vector<int> mMedianHits;
    vector<double> mMeanHits;
    vector<double> mCoverage;
    vector<int> mKmerCounts;
    vector<int> mGenomeReads;
    vector<KCResult> mResults;
    uint32 mNumber;
    uint32 mUniqueHashNum;
    uint32* mHashKCH;
    //KCHit* mKCHits;
    vector<KCHit> mKCHits;
    string mFilename;
    gzFile mZipFile;
    ifstream mFile;
    bool mZipped;
    int mIdBits;
    uint32 mIdMask;
    uint32 mCountMax;
    bool mStatDone;
    uint32 mUniqueNumber;
    //unordered_map<uint64,  KCHit> map_kmer2kch;
    vector<unordered_map<uint64_t, KCHit> > mVec_kh2KCHit;

    mutex mLock;
    data_queue_t mdq;

//change
public:
    unsigned atomic_lock_size;
    unsigned atomic_lock_size_mask;
    //std::atomic_flag *atomic_lock;
    padding_atomic_lock *atomic_lock;
    //vector<KCHit> KCHitmap;
    KCHit *KCHitmap;
    //atomic_ulong mycollect_num;
//change
};


#endif
