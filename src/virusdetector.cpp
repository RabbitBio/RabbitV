#include "virusdetector.h"
#include <atomic>

VirusDetector::VirusDetector(Options* opt){
    mOptions = opt;

    mKmer = NULL;
    if(!mOptions->kmerFile.empty())
        mKmer = new Kmer(mOptions->kmerFile, opt);

    mKmerCollection = NULL;
    if(!mOptions->kmerCollectionFile.empty())
        mKmerCollection = new KmerCollection(mOptions->kmerCollectionFile, opt);

    // no KMER file, the kmerKeyLen is not intialized
    if(mOptions->kmerKeyLen == 0)
        mOptions->kmerKeyLen = 25;
    mGenomes = NULL;
    if(!mOptions->genomeFile.empty())
        mGenomes = new Genomes(mOptions->genomeFile, opt);
    mHits = 0;
    keymask =  (1ul << (2 * mOptions->kmerKeyLen)) - 1;

}

VirusDetector::~VirusDetector(){
    if(mKmer) {
        delete mKmer;
        mKmer = NULL;
    }
    if(mKmerCollection) {
        delete mKmerCollection;
        mKmerCollection = NULL;
    }
    if(mGenomes) {
        delete mGenomes;
        mGenomes = NULL;
    }
}

void VirusDetector::report() {
    if(mKmer) {
        cerr << "Coverage for target unique KMER file:"<<endl;
        mKmer->report();
    }
    if(mKmerCollection) {
        cerr << endl << "Detection result for provided KMER collection:"<<endl;
        mKmerCollection->report();
    }
    if(mGenomes) {
        //mGenomes->report();
    }
}
std::mutex cerr_mutex;
bool VirusDetector::detect(Read* r) {
    if(r->length() >= mOptions->longReadThreshold) {
        // long reads, split it
        vector<Read*> reads = r->split(mOptions->segmentLength);
        bool detected = false;
        for(int i=0; i<reads.size(); i++) {
            // recursive
            detected |= detect(reads[i]);
            delete reads[i];
            reads[i] = NULL;
        }
        return detected;
    }
    string& seq = r->mSeq.mStr;
    //Sequence rSequence = ~(r->mSeq);
    //string& rseq = rSequence.mStr;
    uint32 hitGenomeID = 0;
    bool res = scan(seq, hitGenomeID);
    if( hitGenomeID > 0 ){
      std::lock_guard<std::mutex> lk(cerr_mutex);
      cerr << r->mName << " " << hitGenomeID << " " 
          << mKmerCollection->mNames[hitGenomeID] << endl;
    }

    //return scan(seq); //| scan(rseq);
    return res;
}

inline uint64 max_gid(unordered_map<uint64, int> &gid2count){
  int max = 0;
  uint64 id = 0;
  for (auto x: gid2count){
    //cout << x.first << " - " << x.second << endl;
    if(x.second > max){
      max = x.second;
      id = x.first;
    }
  }
  return id;
}

bool VirusDetector::scan(string& seq, uint32 &hitGenomeID) {
    int hitCount = 0;

    int keylen = mOptions->kmerKeyLen;
    //int blankBits = 64 - 2*keylen;
    unordered_map<uint64, int> gid2count;

    bool onlyHitOneGenome = true;
    uint32 lastGenomeID = 0;

    if(seq.length() < keylen)
        return false;

    bool valid = true;
    bool needAlignment = false;

    uint32 start = 0;
    const uint8_t mask = 0x06;
    uint64 key = Kmer::seq2uint64(seq, start, keylen-1, valid);
    while(valid == false) {
        start++;
        key = Kmer::seq2uint64(seq, start, keylen-1, valid);
        // reach the tail
        if(start >= seq.length() - keylen)
            return false;
    }
    //chang
    int tmplen = seq.length();
    uint32 pos = start + keylen - 1;
    for(; pos < tmplen;) 
    {
        key <<= 2;
        if(seq[pos] == 'N')
        {
            pos++;
            key = 0;
            if(Kmer::seq2uint64(seq, pos, keylen, key, keymask))
                return false;
        }
        else
        {
            uint8_t meri = seq[pos];
            meri &= mask;
            meri >>= 1;
            key |= (uint64_t)meri;
            ++pos;
        }

        key &= keymask;

        if(!needAlignment && mGenomes && mGenomes->hasKey(key))
        {
            needAlignment = true;
        }

        if(mKmer) {
            bool hit = mKmer->add(key);
            if(hit)
                hitCount++;
        }

        if(mKmerCollection) {
          uint32 gid;
          if (mKmerCollection->readBin)
          {
            gid = mKmerCollection->add_bin(key);
            if (gid > 0) gid2count[gid]++;
          }
          else
          {
            gid = mKmerCollection->add(key);
            gid = gid + 1;
          }
          if (gid > 0)
          {
            if (lastGenomeID != 0 && gid != lastGenomeID)
              onlyHitOneGenome = false;
            lastGenomeID = gid;
          }
        }
    }

    if(mKmerCollection && onlyHitOneGenome && lastGenomeID>0){
        mKmerCollection->addGenomeRead(lastGenomeID);
        //hitGenomeID = lastGenomeID;
    }

    hitGenomeID = max_gid(gid2count);
    bool wellMapped = false;
    //if(needAlignment && mGenomes)
    //    wellMapped = mGenomes->align(seq);

    return hitCount>0 || wellMapped;
}
