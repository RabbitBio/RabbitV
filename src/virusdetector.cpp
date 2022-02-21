#include "virusdetector.h"

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

    //because it's reverse complement is included in the unique kmer file or unique kmer collection file,
    //there is no need to detect recerse comolement here.
    return scan(seq); //| scan(rseq);
}

bool VirusDetector::scan(string& seq) {
    int hitCount = 0;

    int keylen = mOptions->kmerKeyLen;
    //int blankBits = 64 - 2*keylen;

    bool onlyHitOneGenome = true;
    uint32 lastGenomeID = 0;

    if(seq.length() < keylen)
        return false;

    bool valid = true;
    bool needAlignment = false;

    uint32 start = 0;
    const uint8_t mask = 0x06;
    //map k-mer to 64-bit int.
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
    //map k-mer to 64-bit int.
    //use special encoding, use bit operations insted of branch operations to imporve running speed.
    // A : 0
    // C : 1
    // T : 2
    // G : 3
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

    if(mKmerCollection && onlyHitOneGenome && lastGenomeID>0)
        mKmerCollection->addGenomeRead(lastGenomeID);

    bool wellMapped = false;
    //if(needAlignment && mGenomes)
    //    wellMapped = mGenomes->align(seq);

    return hitCount>0 || wellMapped;
}
