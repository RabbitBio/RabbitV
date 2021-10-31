#include "kmercollection.h"
#include "util.h"
#include <sstream>
#include <string.h>
#include <memory.h>
#include "kmer.h"
#include <sys/time.h>

inline double get_time(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

const long HASH_LENGTH = (1L<<30);

KmerCollection::KmerCollection(string filename, Options* opt)
{
    double s = get_time();
    mOptions = opt;
    mHashKCH = new uint32[HASH_LENGTH];
    memset(mHashKCH, 0, sizeof(uint32)*HASH_LENGTH);
    mFilename = filename;
    mNumber = 0;
    mIdBits = 0;
    mIdMask = 0;
    mCountMax = 0;
    mStatDone = false;
    mUniqueHashNum = 0;
    //mKCHits = NULL;
    init();
    double e = get_time();
    cerr << "init time : " << e - s << endl;
}

KmerCollection::~KmerCollection()
{
    if(mHashKCH) {
        delete mHashKCH;
        mHashKCH = NULL;
    }

    //if(mKCHits) {
    //    delete mKCHits;
    //    mKCHits = NULL;
    //}

    if (mZipped){
        if (mZipFile){
            gzclose(mZipFile);
            mZipFile = NULL;
        }
    }
    else {
        if (mFile.is_open()){
            mFile.close();
        }
    }
}

bool KmerCollection::getLine(char* line, int maxLine){
    bool status = true;
    if(mZipped)
        status = gzgets(mZipFile, line, maxLine);
    else {
        mFile.getline(line, maxLine);
        status = !mFile.fail();
    }

    // trim \n, \r or \r\n in the tail
    int readed = strlen(line);
    if(readed >=2 ){
        if(line[readed-1] == '\n' || line[readed-1] == '\r'){
            line[readed-1] = '\0';
            if(line[readed-2] == '\r')
                line[readed-2] = '\0';
        }
    }

    return status;
}

bool descComp (int i,int j) { return (i>j); }

bool KCResultComp (KCResult i, KCResult j) { 
    if (i.mCoverage == j.mCoverage)
        return i.mMedianHit > j.mMedianHit;
    else
        return i.mCoverage > j.mCoverage;
}

void KmerCollection::stat(){
    vector<vector<int>> kmerHits(mNumber);
    /*
    for(int i=0; i<mUniqueHashNum; i++) {
        KCHit kch = mKCHits[i];

        if(kch.mHit>0) {
            mHits[kch.mID]+=kch.mHit;
            kmerHits[kch.mID].push_back(kch.mHit);
        }
    }
    */
   for(auto& map_kmer2kch: mVec_kh2KCHit){
    for(auto& k2kch: map_kmer2kch){
      KCHit&  kch = k2kch.second;
      if(kch.mHit > 0){
        mHits[kch.mID] += kch.mHit;
        kmerHits[kch.mID].push_back(kch.mHit);
      }
    }
    }

    for(int id=0; id<mNumber; id++){
        if(mKmerCounts[id] ==  0) {
            mMedianHits[id]=0;
            mMeanHits[id]=0.0;
            mCoverage[id]=0.0;
        } else{
            int medianPos = (mKmerCounts[id]+1)/2;
            if(medianPos >= kmerHits[id].size())
                mMedianHits[id] = 0;
            else {
                nth_element(kmerHits[id].begin(), kmerHits[id].begin()+medianPos, kmerHits[id].end(), descComp);
                mMedianHits[id] = kmerHits[id][medianPos];
            }
            mMeanHits[id] = (double)mHits[id]/(double)mKmerCounts[id];
            mCoverage[id] = (double)kmerHits[id].size()/(double)mKmerCounts[id];
        }
    }

    for(int id=0; id<mNumber; id++){
        if(mCoverage[id] > mOptions->kcCoverageThreshold && mKmerCounts[id] > 10) {
            KCResult kcr;
            kcr.mName = mNames[id];
            kcr.mHit = mHits[id];
            kcr.mCoverage = mCoverage[id];
            kcr.mMedianHit = mMedianHits[id];
            kcr.mMeanHit = mMeanHits[id];
            kcr.mKmerCount = mKmerCounts[id];
            kcr.mUniqueReads = mGenomeReads[id];

            mResults.push_back(kcr);
        }
    }

    sort(mResults.begin(),mResults.end(),KCResultComp);

    mStatDone = true;
}

uint32 KmerCollection::add(uint64 kmer64) {
    uint64 kmerhash = makeHash(kmer64);
    uint32 index = mHashKCH[kmerhash];
    if(index != 0 && index != COLLISION_FLAG) {
        if(mKCHits[index - 1].mKey64 == kmer64) {
            mKCHits[index - 1].mHit++;
            return mKCHits[index - 1].mID + 1;
        } else
            return 0;
    } else
        return 0;
}
uint32 KmerCollection::add_bin(uint64 kmer64) {
    uint64 kmerhash = makeHash(kmer64);
    for(auto& map_kmer2kch : mVec_kh2KCHit){
      auto itr = map_kmer2kch.find(kmerhash);
      if (itr != map_kmer2kch.end())
      {
        if (itr->second.mKey64 = kmer64)
        {
          itr->second.mHit++;
          return itr->second.mID;
        }
      }
    }
    return 0;
}

void KmerCollection::addGenomeRead(uint32 genomeID) {
    if(genomeID-1 < mGenomeReads.size())
        mGenomeReads[genomeID-1]++;
}

void KmerCollection::init()
{
    if(mOptions->verbose)
        loginfo("Initializing k-mer collection: " + mFilename + "\n");
    if (ends_with(mFilename, ".fasta.gz") || ends_with(mFilename, ".fa.gz")){
        mZipFile = gzopen(mFilename.c_str(), "r");
        mZipped = true;
    }
    else if(ends_with(mFilename, ".fasta") || ends_with(mFilename, ".fa")){
        mFile.open(mFilename.c_str(), ifstream::in);
        mZipped = false;
    }else if(ends_with(mFilename, ".binary") || ends_with(mFilename, ".bin")){
        mFile.open(mFilename.c_str(), ifstream::in | ifstream::binary);
        mZipped = false;
        //readAllBin();
        readBin = true;
        mul_thread_init();
        return;
    } else {
        error_exit("Not a FASTA file: " + mFilename);
    }

    //unordered_map<uint64, KCHit> hashKmerMap;

    const int maxLine = 1000;
    char line[maxLine];
    if (mZipped){
        if (mZipFile == NULL)
            return ;
    }

    vector<vector<uint64>> allKmer64;
    int unique = 0;
    int total = 0;
    bool initialized = false;
    while(true) {
        if(eof())
            break;
        getLine(line, maxLine);
        string linestr(line);
        if(linestr.empty() || line[0]=='#') {
            continue;
        }
        if(line[0]=='>') {
            if(total > 0) {
                //cerr << unique << "/" << total << endl;
                mKmerCounts.push_back(unique);
                total=0;
                unique=0;
            }
            mNames.push_back(linestr.substr(1, linestr.length() - 1));
            mHits.push_back(0);
            allKmer64.push_back(vector<uint64>());
            mMeanHits.push_back(0.0);
            mCoverage.push_back(0.0);
            mMedianHits.push_back(0);
            mGenomeReads.push_back(0);
            mNumber++;
            //cerr<<mNumber<<": " << linestr << endl;
            continue;
        }

        string& seq = linestr;
        if(!initialized) {
            initialized = true;
            if(mOptions->kmerKeyLen == 0)
                mOptions->kmerKeyLen = seq.length();

            if(mOptions->kmerKeyLen > 32)
                error_exit("k-mer key length cannot be >32: " + seq);
        }

        if(seq.length() != mOptions->kmerKeyLen) {
            cerr << "k-mer length must be " << mOptions->kmerKeyLen << ", skipped " << seq << endl;
            continue;
        }

        bool valid = true;
        uint64 kmer64 = Kmer::seq2uint64(seq, 0, seq.length(), valid);
        if(valid) {
            total++;
            uint64 kmerhash = makeHash(kmer64);
            //unordered_map<uint64, KCHit>::iterator iter = hashKmerMap.find(kmerhash);
            if(mHashKCH[kmerhash] == 0) {
                unique++;
                mHashKCH[kmerhash] = mNumber;
                allKmer64[mNumber-1].push_back(kmer64);
            } else if(mHashKCH[kmerhash] != COLLISION_FLAG) { // we use the mHits as a flag
                if(mHashKCH[kmerhash] == mNumber)
                    unique--;
                else {
                    int collID = mHashKCH[kmerhash]-1;
                    mKmerCounts[collID]--;
                }
                mHashKCH[kmerhash] = COLLISION_FLAG;
            }
        }
    }
    mKmerCounts.push_back(unique);

    mUniqueHashNum = 0;
    for(int i=0; i<mNumber; i++)
        mUniqueHashNum += mKmerCounts[i];

    //mKCHits = new KCHit[mUniqueHashNum];
    //memset(mKCHits, 0, sizeof(KCHit)*mUniqueHashNum);
    mKCHits = vector<KCHit>(mUniqueHashNum);

    uint32 cur = 0;
    for(int i=0; i<mNumber; i++) {
        for(int j=0; j<allKmer64[i].size(); j++) {
            uint64 kmer64 = allKmer64[i][j];
            uint64 kmerhash = makeHash(kmer64);
            uint32 index = mHashKCH[kmerhash];
            // means unique
            if(index == i+1) {
                if(cur >= mUniqueHashNum)
                    error_exit("Uninque number incorrectly calculated in k-mer collection initialization.");
                mHashKCH[kmerhash] = cur + 1; // 0 means no hit
                mKCHits[cur].mID = i;
                mKCHits[cur].mHit = 0;
                mKCHits[cur].mKey64 = kmer64;
                cur++;
            }
        }
    }

    //makeBitAndMask();
}

bool KmerCollection::eof() {
    if (mZipped) {
        return gzeof(mZipFile);
    } else {
        return mFile.eof();
    }
}

uint64 KmerCollection::makeHash(uint64 key) {
    return (1713137323 * key + (key>>12)*7341234131 + (key>>24)*371371377) & (HASH_LENGTH-1);
}

void KmerCollection::report() {
    if(!mStatDone)
        stat();
    int highConfidenceNum = 0;
    for(int i=0; i<mResults.size(); i++) {
        KCResult kcr = mResults[i];
        if(!isHighConfidence(kcr))
            continue;
        highConfidenceNum++ ;
        cerr << highConfidenceNum << ", ";
        cerr << kcr.mName << ",";
        cerr << "coverage:" << kcr.mCoverage;
        cerr << ",kmer_count:" << kcr.mKmerCount;
        cerr << ",median_depth:" << kcr.mMedianHit;
        cerr << ",mean_depth:" << kcr.mMeanHit;
        cerr << ",unique_reads:" << kcr.mUniqueReads;
        cerr <<  endl;
    }
    if(highConfidenceNum == 0)
        cerr << "No high confidence k-mer coverage found." << endl;
}

void KmerCollection::reportJSON(ofstream& ofs) {
    if(!mStatDone)
        stat();

    ofs << "\t" << "\"kmer_collection_scan_result\": {" << endl;

    int first = true;
    for(int i=0; i<mResults.size(); i++) {
        KCResult kcr = mResults[i];
        if(first) {
            first = false;
        } else 
            ofs << "," << endl;
        string name = replace(kcr.mName, "\"", "'");
        ofs << "\t\t\"" << name << "\":{";
        ofs << "\"coverage\":" << kcr.mCoverage;
        ofs << ",\"kmer_count\":" << kcr.mKmerCount;
        ofs << ",\"kmer_hits\":" << kcr.mHit;
        ofs << ",\"median_depth\":" << kcr.mMedianHit;
        ofs << ",\"mean_depth\":" << kcr.mMeanHit;
        ofs << ",\"unique_reads\":" << kcr.mUniqueReads;
        ofs << "}";
    }

    ofs << endl << "\t}," << endl;
}

bool KmerCollection::isHighConfidence(KCResult kcr) {
    if(kcr.mCoverage > mOptions->kcCoverageHighConfidence && kcr.mMedianHit > mOptions->kcMedianHitHighConfidence)
        return true;
    else
        return false;
}

void KmerCollection::reportHTML(ofstream& ofs) {
    ofs << "<table class='summary_table' style='width:100%'>\n";
    ofs <<  "<tr style='background:#cccccc'> <td>Genome</td><td>K-mer hits</td><td>Unique reads</td><td>Coverage</td><td>Median depth</td><td>Mean depth</td><td>Remark</td>  </tr>"  << endl;

    int highConfidence = 0;
    for(int i=0; i<mResults.size(); i++) {
        KCResult kcr = mResults[i];

        if(!isHighConfidence(kcr))
            continue;

        highConfidence++;

        string remark;
        if(kcr.mName.find("phi-X174")!=string::npos)
            remark = "Illumina PhiX control library";
        else if(kcr.mName ==  "NC_045512.2 Wuhan seafood market pneumonia virus isolate Wuhan-Hu-1, complete genome")
            remark = "SARS-CoV-2";
        ofs << "<tr>";
        ofs << "<td width=52%>" << kcr.mName << " (" << kcr.mKmerCount << " k-mer keys)</td>";
        ofs << "<td width=8%>" << kcr.mHit << "</td>";
        ofs << "<td width=8%>" << kcr.mUniqueReads << "</td>";
        ofs << "<td width=10%>" << kcr.mCoverage * 100 << "%</td>";
        ofs << "<td width=8%>" << kcr.mMedianHit << "</td>";
        ofs << "<td width=8%>" << kcr.mMeanHit << "</td>";
        ofs << "<td width=6%>" << remark << "</td>";
        ofs << "</tr>" <<  endl;
    }
    if(highConfidence == 0)
        ofs << "<tr> <td colspan=6 style='text-align:center;'>No high confidence k-mer coverage found. </td></tr>" << endl;
    ofs << "</table>\n";

    if(highConfidence != mResults.size())  {
        ofs << "<div class='subsection_title' style='font-size:12px;font-weight:normal;color:#223399;' onclick=showOrHide('low_confidence_kcr')>+ Show " << mResults.size() - highConfidence;
        ofs << " more with low confidence (coverage <= " << mOptions->kcCoverageHighConfidence *  100;
        ofs << "% or median depth <= " << mOptions->kcMedianHitHighConfidence;
        ofs << ") ▼ </div>\n";
        ofs << "<table id='low_confidence_kcr' style='display:none;width:100%;' class='summary_table'>\n";
        ofs <<  "<tr style='background:#cccccc'> <td>Genome</td><td>K-mer hits</td><td>Unique reads</td><td>Coverage</td><td>Median depth</td><td>Mean depth</td><td>Remark</td>  </tr>"  << endl;

        int highConfidence = 0;
        for(int i=0; i<mResults.size(); i++) {
            KCResult kcr = mResults[i];

            if(isHighConfidence(kcr))
                continue;

            string remark;
            if(kcr.mName.find("phi-X174")!=string::npos)
                remark = "Illumina PhiX control library";
            else if(kcr.mName ==  "NC_045512.2 Wuhan seafood market pneumonia virus isolate Wuhan-Hu-1, complete genome")
                remark = "SARS-CoV-2";
            ofs << "<tr>";
            ofs << "<td width=52%>" << kcr.mName << " (" << kcr.mKmerCount << " k-mer keys)</td>";
            ofs << "<td width=8%>" << kcr.mHit << "</td>";
            ofs << "<td width=8%>" << kcr.mUniqueReads << "</td>";
            ofs << "<td width=10%>" << kcr.mCoverage * 100 << "%</td>";
            ofs << "<td width=8%>" << kcr.mMedianHit << "</td>";
            ofs << "<td width=8%>" << kcr.mMeanHit << "</td>";
            ofs << "<td width=6%>" << remark << "</td>";
            ofs << "</tr>" <<  endl;
        }
        ofs << "</table>\n";
    }
}

void KmerCollection::makeBitAndMask() {
    int bits=1;
    uint32 mask =0x01;

    while(true){
        if(bits == 32) {
            error_exit("Too many contigs in: " + mFilename);
        }
        if(mask > mNumber)
            break;
        bits++;
        mask= (mask << 1) + 1;
    }
    mIdBits = bits;
    mIdMask = mask;
    mCountMax  = (0xFFFFFFFFFFF >> mIdBits);
}

uint32 KmerCollection::packIdCount(uint32 id, uint32 count) {
    uint32 data = 0;
    data |= (count << mIdBits); // high bits: count
    data |= id; // low bits: id
    return data;
}

void KmerCollection::unpackIdCount(uint32 data, uint32& id, uint32& count) {
    count = data >> mIdBits;
    id = data & mIdMask;
}

void KmerCollection::readAllBin() {
  //unordered_map<uint32, uint64> id2enar;//id2enar -> id to en array number mapper
  vector<uint64> id2enar;
  char fname[256];
  int unique = 0;
  int total = 0;
  vector<vector<uint64_t>> allKmer64;
  string current_ref;
  if (mOptions->kmerKeyLen == 0 || mOptions->kmerKeyLen > 32)
    cerr << "plase specify a <= 32 kmer length (default = 25 )"  << endl;
  uint64 coll_num = 0;
  while(!mFile.eof()){
    //-- reading name
    mFile.getline(fname, 256);
    current_ref = string(fname);
    if(current_ref == ""){
      continue;
    }
    if(total > 0){
      mKmerCounts.push_back(unique);
      total = 0; unique = 0;
    }
    mNames.push_back(current_ref);
    mHits.push_back(0);
    allKmer64.push_back(vector<uint64>());
    mMeanHits.push_back(0.0);
    mCoverage.push_back(0.0);
    mMedianHits.push_back(0);
    mGenomeReads.push_back(0);
    id2enar.push_back(0);
    mNumber++;

    //-- read content
    uint64_t ksize;
    mFile.read(reinterpret_cast<char*>(&ksize), 8);
    cerr << mNumber << ": " << current_ref << " kmer size: " << ksize << endl;
    uint64_t *kmers = new uint64_t[ksize];
    mFile.read(reinterpret_cast<char*>(kmers), ksize * 8);
    for(int i = 0; i < ksize; ++i){
      uint64_t kmer64 = kmers[i];
      total++;
      uint64 kmerhash = makeHash(kmer64);
      //unordered_map<uint64, KCHit>::iterator iter = hashKmerMap.find(kmerhash);
      if (mHashKCH[kmerhash] == 0)
      {
        unique++;
        //mHashKCH[kmerhash] = mNumber;
        //allKmer64[mNumber - 1].push_back(kmer64);
        //-----------------------
        //if (curr >= unique)
        //  error_exit("Uninque number incorrectly calculated in k-mer collection initialization.");
        //mHashKCH[kmerhash] = mNumber; // 0 means no hit
        //mKCHits[cur].mID = i;
        //mKCHits[cur].mHit = 0;
        //mKCHits[cur].mKey64 = kmer64;
        mKCHits.push_back(KCHit{kmer64, mNumber-1, 0});
        mUniqueHashNum++;
        mHashKCH[kmerhash] = mUniqueHashNum; // 0 means no hit
        id2enar[mNumber-1]++;
        //-----------------------
      }
      else{
        coll_num++;
        KCHit& coll_kch = mKCHits[mHashKCH[kmerhash]];
        const uint32 coll_id = coll_kch.mID;
        const uint32 my_id = mNumber-1;
        if(id2enar[my_id] < id2enar[coll_id]+1){ 
          coll_kch = KCHit{kmer64, mNumber-1, 0};
          id2enar[my_id]++;
          id2enar[coll_id]--;
        } 
      }
      //else if (mHashKCH[kmerhash] != COLLISION_FLAG)
      //{ // we use the mHits as a flag
      //  if (mHashKCH[kmerhash] == mNumber){
      //    unique--;
      //    mUniqueHashNum--;
      //  }else{
      //    int collID = mHashKCH[kmerhash] - 1;
      //    //cout << "collID " << collID << endl;
      //    //cout << "count: " << mKmerCounts[collID] << endl;
      //    //mKmerCounts[collID]--;
      //  }
      //  mHashKCH[kmerhash] = COLLISION_FLAG;
      //}
    }
    delete [] kmers;
  }
  cout << "write over! " << "conflict number: " << coll_num 
      << " and unique number: " << mUniqueHashNum << endl; 
  /*
  mKmerCounts.push_back(unique);
  mKmerCounts.push_back(unique);
  mUniqueHashNum = 0;
  for (int i = 0; i < mNumber; i++)
    mUniqueHashNum += mKmerCounts[i];

  mKCHits = new KCHit[mUniqueHashNum];
  memset(mKCHits, 0, sizeof(KCHit) * mUniqueHashNum);

  uint32 cur = 0;
  for (int i = 0; i < mNumber; i++)
  {
    for (int j = 0; j < allKmer64[i].size(); j++)
    {
      uint64 kmer64 = allKmer64[i][j];
      uint64 kmerhash = makeHash(kmer64);
      uint32 index = mHashKCH[kmerhash];
      // means unique
      if (index == i + 1)
      {
        if (cur >= mUniqueHashNum)
          error_exit("Uninque number incorrectly calculated in k-mer collection initialization.");
        mHashKCH[kmerhash] = cur + 1; // 0 means no hit
        mKCHits[cur].mID = i;
        mKCHits[cur].mHit = 0;
        mKCHits[cur].mKey64 = kmer64;
        cur++;
      }
    }
  }
  */
}

void KmerCollection::mul_thread_init(){
  auto lambda_rf = [&]() -> void
  {
    vector<uint64> id2enar;
    char fname[256];
    int unique = 0;
    int total = 0;
    string current_ref;
    if (mOptions->kmerKeyLen == 0 || mOptions->kmerKeyLen > 32)
      cerr << "plase specify a <= 32 kmer length (default = 25 )" << endl;
    uint64 coll_num = 0;
    while (!mFile.eof())
    {
      //-- reading name
      mFile.getline(fname, 256);
      current_ref = string(fname);
      if (current_ref == "")
      {
        continue;
      }
      //if (total > 0)
      //{
      //  mKmerCounts.push_back(unique);
      //  total = 0;
      //  unique = 0;
      //}
      mKmerCounts.push_back(0);
      mNames.push_back(current_ref);
      mHits.push_back(0);
      mMeanHits.push_back(0.0);
      mCoverage.push_back(0.0);
      mMedianHits.push_back(0);
      mGenomeReads.push_back(0);
      id2enar.push_back(0);
      mNumber++;

      //-- read content
      uint64_t ksize;
      mFile.read(reinterpret_cast<char *>(&ksize), 8);
      //{
      //  lock_guard<mutex> lk(mLock);
      //  cerr << static_cast<void *>(this) << " " << mNumber << ": " << current_ref << " kmer size: " << ksize << endl;
      //}
      uint64_t *kmers = new uint64_t[ksize];
      mFile.read(reinterpret_cast<char *>(kmers), ksize * 8);
      auto tmp = new std::tuple<uint32, uint64_t*, uint64_t>{mNumber - 1, kmers, ksize};
      mdq.Push(mNumber, tmp);
    }
    //mFinished = true;
    mdq.SetCompleted();
  };
  auto lambda_wf = [&](unordered_map<uint64_t, KCHit> &map_kmer2kch) -> void
  {
    //unordered_map<uint64_t, KC_t> map_kmer2kch;
    std::tuple<uint32, uint64_t *, uint64_t> *data;// = mdq.front();
    rabbit::int64 uid;
    while (mdq.Pop(uid, data))
    {
      uint32 ref_id;
      uint64_t unique = 0;
      uint64_t *kmers;
      uint64_t ksize;
      ref_id = std::get<0>(*data);
      kmers = std::get<1>(*data);
      ksize = std::get<2>(*data);
      //cout << "info : " << static_cast<void *>(this) << " " << ref_id << " " << kmers << " " << ksize << endl;
      for (int i = 0; i < ksize; ++i)
      {
        uint64_t kmer64 = kmers[i];
        uint64 kmerhash = makeHash(kmer64);
        if (!map_kmer2kch.count(kmerhash))
        {
          map_kmer2kch[kmerhash] = KCHit{kmer64, ref_id, 0};
          unique++;
          //map_kmer2kch[kmer64] = KC_t{ref_id, 0};
        }
        //else
        //{
        //}
      }
      mKmerCounts[ref_id] = unique;
      delete[] kmers;
    }
  };
  std::thread rt(lambda_rf); //read function -> 1
  vector<std::thread> wt;    //worker function -> nth
  mVec_kh2KCHit = vector<unordered_map<uint64_t, KCHit> >(mOptions->thread);
  double t1 = get_time();
  for(int i = 0; i < mOptions->thread; i++){
    wt.push_back(std::thread(std::bind(lambda_wf, std::ref(mVec_kh2KCHit[i]))));
  }
  rt.join();
  for(auto &t : wt)
    t.join();
  double t2 = get_time();

  //for(int i = 1; i < maps.size(); i++){
  //  maps[0].insert(maps[i].begin(), maps[i].end());
  //  unordered_map<uint64_t, KCHit>().swap(maps[i]);
  //}

  double t3 = get_time();
  //this->map_kmer2kch = std::move(maps[0]);
  cout << "first step time: " << t2 - t1 
        << " second step time: " << t3 - t2
        << " total time: " << t3 - t1 << endl;
  //process mKmerCount

}