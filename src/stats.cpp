#include "stats.h"
#include <memory.h>
#include <sstream>
#include <immintrin.h>
#include "util.h"

#define KMER_LEN 5

Stats::Stats(Options* opt, bool isRead2, int guessedCycles, int bufferMargin){
    mOptions = opt;
    mIsRead2 = isRead2;
    mReads = 0;
    mLengthSum = 0;

    mEvaluatedSeqLen = mOptions->seqLen1;
    if(mIsRead2)
        mEvaluatedSeqLen = mOptions->seqLen2;

    if(guessedCycles == 0) {
        guessedCycles = mEvaluatedSeqLen;
    }

    mCycles = guessedCycles;
    mBases = 0;
    mQ20Total = 0;
    mQ30Total = 0;
    summarized = false;
    mKmerMin = 0;
    mKmerMax = 0;

    // extend the buffer to make sure it's long enough
    mBufLen = guessedCycles + bufferMargin;

    for(int i=0; i<8; i++){
        mQ20Bases[i] = 0;
        mQ30Bases[i] = 0;
        mBaseContents[i] = 0;

        mCycleQ30Bases[i] = new long[mBufLen];
        memset(mCycleQ30Bases[i], 0, sizeof(long) * mBufLen);

        mCycleQ20Bases[i] = new long[mBufLen];
        memset(mCycleQ20Bases[i], 0, sizeof(long) * mBufLen);

        mCycleBaseContents[i] = new long[mBufLen];
        memset(mCycleBaseContents[i], 0, sizeof(long) * mBufLen);

        mCycleBaseQual[i] = new long[mBufLen];
        memset(mCycleBaseQual[i], 0, sizeof(long) * mBufLen);
    }
    mCycleTotalBase = new long[mBufLen];
    memset(mCycleTotalBase, 0, sizeof(long)*mBufLen);

    mCycleTotalQual = new long[mBufLen];
    memset(mCycleTotalQual, 0, sizeof(long)*mBufLen);

    //临时数组
    mCycleQ30BasesTemp = new int[mBufLen*8];
    memset(mCycleQ30BasesTemp, 0, sizeof(int)*mBufLen*8);

    mCycleQ20BasesTemp = new int[mBufLen*8];
    memset(mCycleQ20BasesTemp, 0, sizeof(int)*mBufLen*8);

    mCycleBaseContentsTemp = new int[mBufLen*8];
    memset(mCycleBaseContentsTemp, 0, sizeof(int)*mBufLen*8);

    mCycleBaseQualTemp = new int[mBufLen*8];
    memset(mCycleBaseQualTemp, 0, sizeof(int)*mBufLen*8);
    //将临时数组剩于空间序列条数初始化为maxCount
    tempEmpty = maxCount;

    /**没有必要有mkmer了
    mKmerBufLen = 2<<(KMER_LEN * 2);
    mKmer = new long[mKmerBufLen];
    memset(mKmer, 0, sizeof(long)*mKmerBufLen);
    **/
}

void Stats::extendBuffer(int newBufLen){
    if(newBufLen <= mBufLen)
        return ;

    long* newBuf = NULL;

    for(int i=0; i<8; i++){
        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleQ30Bases[i], sizeof(long) * mBufLen);
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleQ20Bases[i], sizeof(long) * mBufLen);
        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleBaseContents[i], sizeof(long) * mBufLen);
        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleBaseQual[i], sizeof(long) * mBufLen);
        delete mCycleBaseQual[i];
        mCycleBaseQual[i] = newBuf;
    }
    newBuf = new long[newBufLen];
    memset(newBuf, 0, sizeof(long)*newBufLen);
    memcpy(newBuf, mCycleTotalBase, sizeof(long)*mBufLen);
    delete mCycleTotalBase;
    mCycleTotalBase = newBuf;

    newBuf = new long[newBufLen];
    memset(newBuf, 0, sizeof(long)*newBufLen);
    memcpy(newBuf, mCycleTotalQual, sizeof(long)*mBufLen);
    delete mCycleTotalQual;
    mCycleTotalQual = newBuf;

    //临时数组的扩展
    int *newTempBuff = NULL;

    newTempBuff = new int[newBufLen*8];
    memset(newTempBuff, 0, sizeof(int)*newBufLen*8);
    memcpy(newTempBuff, mCycleQ30BasesTemp, sizeof(int)*mBufLen*8);
    delete mCycleQ30BasesTemp;
    mCycleQ30BasesTemp = newTempBuff;

    newTempBuff = new int[newBufLen*8];
    memset(newTempBuff, 0, sizeof(int)*newBufLen*8);
    memcpy(newTempBuff, mCycleQ20BasesTemp, sizeof(int)*mBufLen*8);
    delete mCycleQ20BasesTemp;
    mCycleQ20BasesTemp = newTempBuff;

    newTempBuff = new int[newBufLen*8];
    memset(newTempBuff, 0, sizeof(int)*newBufLen*8);
    memcpy(newTempBuff, mCycleBaseContentsTemp, sizeof(int)*mBufLen*8);
    delete mCycleBaseContentsTemp;
    mCycleBaseContentsTemp = newTempBuff;

    newTempBuff = new int[newBufLen*8];
    memset(newTempBuff, 0, sizeof(int)*newBufLen*8);
    memcpy(newTempBuff, mCycleBaseQualTemp, sizeof(int)*mBufLen*8);
    delete mCycleBaseQualTemp;
    mCycleBaseQualTemp = newTempBuff;


    mBufLen = newBufLen;
}

Stats::~Stats() {
    for(int i=0; i<8; i++){
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = NULL;

        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = NULL;

        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = NULL;

        delete mCycleBaseQual[i];
        mCycleBaseQual[i] = NULL;
    }

    delete mCycleTotalBase;
    delete mCycleTotalQual;

    //临时数组的删除
    delete mCycleQ30BasesTemp;
    delete mCycleQ20BasesTemp;
    delete mCycleBaseContentsTemp;
    delete mCycleBaseQualTemp;

    // delete memory of curves
    map<string, double*>::iterator iter;
    for(iter = mQualityCurves.begin(); iter != mQualityCurves.end(); iter++) {
        delete iter->second;
    }
    for(iter = mContentCurves.begin(); iter != mContentCurves.end(); iter++) {
        delete iter->second;
    }
    /**没有mkmer了
    delete mKmer;
    **/
}

void Stats::summarize(bool forced) {
    if(summarized && !forced)
        return;

    moveTempArr2FinalArr();//将临时数组中的值移入最终数组

    countmCycleTotalBase();//根据各分base计算总base
    countmCycleTotalQual();//根据各分qual计算总qual

    // first get the cycle and count total bases
    for(int c=0; c<mBufLen; c++) {
        mBases += mCycleTotalBase[c];
        if (mCycleTotalBase[c] == 0){
            mCycles = c;
            break;
        }
    }
    if(mCycleTotalBase[mBufLen-1]>0)
        mCycles = mBufLen;

    // Q20, Q30, base content
    for(int i=0; i<8; i++) {
        for(int c=0; c<mCycles; c++) {
            mQ20Bases[i] += mCycleQ20Bases[i][c];
            mQ30Bases[i] += mCycleQ30Bases[i][c];
            mBaseContents[i] += mCycleBaseContents[i][c];
        }
        mQ20Total += mQ20Bases[i];
        mQ30Total += mQ30Bases[i];
    }


    // quality curve for mean qual
    double* meanQualCurve = new double[mCycles];
    memset(meanQualCurve, 0, sizeof(double)*mCycles);
    for(int c=0; c<mCycles; c++) {
        meanQualCurve[c] = (double)mCycleTotalQual[c] / (double)mCycleTotalBase[c];
    }
    mQualityCurves["mean"] = meanQualCurve;

    // quality curves and base content curves for different nucleotides
    char alphabets[5] = {'A', 'C', 'T', 'G', 'N'};
    for(int i=0; i<5; i++) {
        char base = alphabets[i];
        // get last 3 bits
        char b = base & 0x07;
        double* qualCurve = new double[mCycles];
        memset(qualCurve, 0, sizeof(double)*mCycles);
        double* contentCurve = new double[mCycles];
        memset(contentCurve, 0, sizeof(double)*mCycles);
        for(int c=0; c<mCycles; c++) {
            if(mCycleBaseContents[b][c] == 0)
                qualCurve[c] = meanQualCurve[c];
            else
                qualCurve[c] = (double)mCycleBaseQual[b][c] / (double)mCycleBaseContents[b][c];
            contentCurve[c] = (double)mCycleBaseContents[b][c] / (double)mCycleTotalBase[c];
        }
        mQualityCurves[string(1, base)] = qualCurve;
        mContentCurves[string(1, base)] = contentCurve;
    }

    // GC content curve
    double* gcContentCurve = new double[mCycles];
    memset(gcContentCurve, 0, sizeof(double)*mCycles);
    char gBase = 'G' & 0x07;
    char cBase = 'C' & 0x07;
    for(int c=0; c<mCycles; c++) {
        gcContentCurve[c] = (double)(mCycleBaseContents[gBase][c] + mCycleBaseContents[cBase][c]) / (double)mCycleTotalBase[c];
    }
    mContentCurves["GC"] = gcContentCurve;

    /**把涉及到mkmer的没必要的部分删去
    mKmerMin = mKmer[0];
    mKmerMax = mKmer[0];
    for(int i=0; i<mKmerBufLen; i++) {
        if(mKmer[i] > mKmerMax)
            mKmerMax = mKmer[i];
        if(mKmer[i] < mKmerMin)
            mKmerMin = mKmer[i];
    }
    **/

    summarized = true;
}

int Stats::getMeanLength() {
    if(mReads == 0)
        return 0.0;
    else
        return mLengthSum/mReads;
}

void Stats::countmCycleTotalBase(){
    //这里将mCycleTotalBase计算出来
    for(int i = 0;i<mBufLen;i++){
        mCycleTotalBase[i] = 0;
        for(int j = 0;j<8;j++){
            mCycleTotalBase[i] += mCycleBaseContents[j][i];
        }
    }
}

void Stats::countmCycleTotalQual(){
    //这里将mCycleTotalQual计算出来
    for(int i = 0;i<mBufLen;i++){
        mCycleTotalQual[i] = 0;
        for(int j = 0;j<8;j++){
            mCycleTotalQual[i] += mCycleBaseQual[j][i];
        }
    }
}

void Stats::moveTempArr2FinalArr(){
    //将临时数组里的东西移到最终数组当中，并且将临时数组归零
    for(int i = 0;i<mBufLen;i++){
        for(int j = 0;j<8;j++){
            mCycleQ30Bases[j][i] += mCycleQ30BasesTemp[i*8+j];
            mCycleQ20Bases[j][i] += mCycleQ20BasesTemp[i*8+j];
            mCycleBaseContents[j][i] += mCycleBaseContentsTemp[i*8+j];
            mCycleBaseQual[j][i] += mCycleBaseQualTemp[i*8+j];

            mCycleQ30BasesTemp[i*8+j] = 0;
            mCycleQ20BasesTemp[i*8+j] = 0;
            mCycleBaseContentsTemp[i*8+j] = 0;
            mCycleBaseQualTemp[i*8+j] = 0;
        }
    }
    tempEmpty = maxCount;
}

#if defined __AVX512F__ && defined __AVX512VL__ && defined __AXV512BITALG__ && defined __MAVX512BW__ && defined __MAXV512VBMI2__

void Stats::statRead(Read* r) 
{
    int len = r->length();
    if(len>268435455){
        //因为mCycleQ30BasesTemp之类的临时数组是buffLen长度的8倍，
        //而要想用int作为下标进行索引，就要求序列的长度要小于(2^31-1)/8=268435455
        cerr<<"序列太长啦，最多只能支持268435455长度的序列"<<endl;
    }

    mLengthSum += len;

    if(mBufLen < len) {
        extendBuffer(max(len + 100, (int)(len * 1.5)));
    }
    const char* seqstr = r->mSeq.mStr.c_str();
    const char* qualstr = r->mQuality.c_str();

    int pallTime = len/16;//向量化最多循环次数
    int leftTime = len - pallTime*16;//剩于的不能向量化的部分

    const int constInt0X07 = (int)(0X07);
    const __m512i const0X07 = _mm512_set1_epi32(constInt0X07);//创建一个0X07常数向量
    const __m512i constAdd1 = _mm512_set1_epi32(1);//创建一个用于+1的常数向量
    const __m512i constSub33 = _mm512_set1_epi32(33);//创建一个用于-33的常数向量
    const __m512i constAdd128 = _mm512_set1_epi32(128);//创建一个用于+128的常数向量
    const __m512i constCharq20 = _mm512_set1_epi32('5');//创建一个q20常数向量
    const __m512i constCharq30 = _mm512_set1_epi32('?');//创建一个q30常数向量

    __m512i storeBaseIndex = _mm512_set_epi32(120, 112, 104, 96, 88, 80, 72, 64, 56, 48, 40, 32, 24, 16, 8, 0);
    
    for(int i = 0;i<pallTime;i++){
        __mmask16 mask = (1ul << 16) - 1;//16位一次全部load进来
        __m128i base128 = _mm_maskz_loadu_epi8(mask, seqstr+i*16);//load到一个128空间，一共有16个char
        //然后让每个char扩展到32位，使用的ZeroExtend32
        __m512i base512 = _mm512_cvtepu8_epi32(base128);

        __m128i qual128 = _mm_maskz_loadu_epi8(mask, qualstr+i*16);
        __m512i qual512 = _mm512_cvtepu8_epi32(qual128);//这个使用的是0扩展的，ZeroExtend32

        __m512i offset = _mm512_and_epi32(base512, const0X07);
        __m512i storeIndex = _mm512_add_epi32(storeBaseIndex, offset);

        __mmask16 q30mask = _mm512_cmp_epi32_mask(constCharq30, qual512, _MM_CMPINT_LE);//得到大于等于q30的mask
        __mmask16 q20mask = _mm512_cmp_epi32_mask(constCharq20, qual512, _MM_CMPINT_LE);//得到大于等于q20的mask

        //把mCycleQ30BasesTemp增加完毕
        __m512i mCycleQ30BasesTempGather = _mm512_i32gather_epi32(storeIndex, mCycleQ30BasesTemp, sizeof(int));
        __m512i mCycleQ30BasesTempStore = _mm512_mask_add_epi32(mCycleQ30BasesTempGather, q30mask, mCycleQ30BasesTempGather, constAdd1);
        _mm512_i32scatter_epi32(mCycleQ30BasesTemp, storeIndex, mCycleQ30BasesTempStore, sizeof(int));

        //把mCycleQ20BasesTemp增加完毕
        __m512i mCycleQ20BasesTempGather = _mm512_i32gather_epi32(storeIndex, mCycleQ20BasesTemp, sizeof(int));
        __m512i mCycleQ20BasesTempStore = _mm512_mask_add_epi32(mCycleQ20BasesTempGather, q20mask, mCycleQ20BasesTempGather, constAdd1);
        _mm512_i32scatter_epi32(mCycleQ20BasesTemp, storeIndex, mCycleQ20BasesTempStore, sizeof(int));

        //把mCycleBaseContentsTemp增加完毕
        __m512i mCycleBaseContentsTempGather = _mm512_i32gather_epi32(storeIndex, mCycleBaseContentsTemp, (int)sizeof(int));
        __m512i mCycleBaseContentsTempStore = _mm512_add_epi32(mCycleBaseContentsTempGather, constAdd1);
        _mm512_i32scatter_epi32(mCycleBaseContentsTemp, storeIndex, mCycleBaseContentsTempStore, (int)sizeof(int));

        //获得需要叠加上去的质量
        __m512i addQual = _mm512_sub_epi32(qual512, constSub33);
        __m512i mCycleBaseQualTempGather = _mm512_i32gather_epi32(storeIndex, mCycleBaseQualTemp, sizeof(int));
        __m512i mCycleBaseQualTempStore = _mm512_add_epi32(mCycleBaseQualTempGather, addQual);
        _mm512_i32scatter_epi32(mCycleBaseQualTemp, storeIndex, mCycleBaseQualTempStore, sizeof(int));

        __m512i temp = _mm512_add_epi32(storeBaseIndex, constAdd128);
        storeBaseIndex = temp;
    }

    //把剩下的给一个个算了
    for(int i = pallTime*16; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        if(qual >= q30) {
            mCycleQ30BasesTemp[8*i+b]++;
            mCycleQ20BasesTemp[8*i+b]++;
        } else if(qual >= q20) {
            mCycleQ20BasesTemp[8*i+b]++;
        }

        mCycleBaseContentsTemp[8*i+b]++;
        mCycleBaseQualTemp[8*i+b] += (qual-33);

        /**
        mCycleTotalBase[i]++;//这个东西可以最后通过mCycleBaseContents累加起来得到，所以没有必要在这里计算
        mCycleTotalQual[i] += (qual-33);//这个东西可以最后通过mCycleBaseQual累加起来得到，所以也没有必要在这里计算
        **/
    }

    tempEmpty--;
    if(tempEmpty<=0){
        //如果剩于容量都用完了，就要把临时数组里面的内容搬运到最终数组当中，并且将临时数组归零
        moveTempArr2FinalArr();
    }

    mReads++;
}

#else
void Stats::statRead(Read* r) {
    int len = r->length();

    mLengthSum += len;

    if(mBufLen < len) {
        extendBuffer(max(len + 100, (int)(len * 1.5)));
    }
    const char* seqstr = r->mSeq.mStr.c_str();
    const char* qualstr = r->mQuality.c_str();

    int kmer = 0;
    bool needFullCompute = true;
    for(int i=0; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        if(qual >= q30) {
            mCycleQ30Bases[b][i]++;
            mCycleQ20Bases[b][i]++;
        } else if(qual >= q20) {
            mCycleQ20Bases[b][i]++;
        }

        mCycleBaseContents[b][i]++;
        mCycleBaseQual[b][i] += (qual-33);

        mCycleTotalBase[i]++;
        mCycleTotalQual[i] += (qual-33);

        if(base == 'N'){
            needFullCompute = true;
            continue;
        }

        // 5 bases required for kmer computing
        if(i<4)
            continue;

        // calc 5 KMER
        // 0x3FC == 0011 1111 1100
        if(!needFullCompute){
            int val = base2val(base);
            if(val < 0){
                needFullCompute = true;
                continue;
            } else {
                kmer = ((kmer<<2) & 0x3FC ) | val;
                //mKmer[kmer]++;
            }
        } else {
            bool valid = true;
            kmer = 0;
            for(int k=0; k<5; k++) {
                int val = base2val(seqstr[i - 4 + k]);
                if(val < 0) {
                    valid = false;
                    break;
                }
                kmer = ((kmer<<2) & 0x3FC ) | val;
            }
            if(!valid) {
                needFullCompute = true;
                continue;
            } else {
                //mKmer[kmer]++;
                needFullCompute = false;
            }
        }

    }

    mReads++;
}

#endif

int Stats::base2val(char base) {
    // (A & 0X07)>>1 = （000）
    // (C & 0X07)>>1 = （001）
    // (G & 0X07)>>1 = （011）
    // (T & 0X07)>>1 = （010）
    if(base=='N'){
        return -1;
    }else{
        return (base & 0X07)>>1;
    }
    // switch(base){
    //     case 'A':
    //         return 0;
    //     case 'T':
    //         return 2;
    //     case 'C':
    //         return 1;
    //     case 'G':
    //         return 3;
    //     default:
    //         return -1;
    // }
}

int Stats::getCycles() {
    if(!summarized)
        summarize();
    return mCycles;
}

long Stats::getReads() {
    if(!summarized)
        summarize();
    return mReads;
}

long Stats::getBases() {
    if(!summarized)
        summarize();
    return mBases;
}

long Stats::getQ20() {
    if(!summarized)
        summarize();
    return mQ20Total;
}

long Stats::getQ30() {
    if(!summarized)
        summarize();
    return mQ30Total;
}

long Stats::getGCNumber() {
    if(!summarized)
        summarize();
    return mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07];
}

void Stats::print() {
    if(!summarized) {
        summarize();
    }
    cerr << "total reads: " << mReads << endl;
    cerr << "total bases: " << mBases << endl;
    cerr << "Q20 bases: " << mQ20Total << "(" << (mQ20Total*100.0)/mBases << "%)" << endl;
    cerr << "Q30 bases: " << mQ30Total << "(" << (mQ30Total*100.0)/mBases << "%)" << endl;
}

void Stats::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"total_reads\": " << mReads << "," << endl;
    ofs << padding << "\t" << "\"total_bases\": " << mBases << "," << endl;
    ofs << padding << "\t" << "\"q20_bases\": " << mQ20Total << "," << endl;
    ofs << padding << "\t" << "\"q30_bases\": " << mQ30Total << "," << endl;
    ofs << padding << "\t" << "\"total_cycles\": " << mCycles << "," << endl;

    // quality curves
    string qualNames[5] = {"A", "T", "C", "G", "mean"};
    ofs << padding << "\t" << "\"quality_curves\": {" << endl;
    for(int i=0 ;i<5; i++) {
        string name=qualNames[i];
        double* curve = mQualityCurves[name];
        ofs << padding << "\t\t" << "\"" << name << "\":[";
        for(int c = 0; c<mCycles; c++) {
            ofs << curve[c];
            // not the end
            if(c != mCycles - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(i != 5-1)
            ofs << ",";
        ofs << endl; 
    }
    ofs << padding << "\t" << "}," << endl;

    // content curves
    string contentNames[6] = {"A", "T", "C", "G", "N", "GC"};
    ofs << padding << "\t" << "\"content_curves\": {" << endl;
    for(int i=0 ;i<6; i++) {
        string name=contentNames[i];
        double* curve = mContentCurves[name];
        ofs << padding << "\t\t" << "\"" << name << "\":[";
        for(int c = 0; c<mCycles; c++) {
            ofs << curve[c];
            // not the end
            if(c != mCycles - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(i != 6-1)
            ofs << ",";
        ofs << endl; 
    }
    ofs << padding << "\t" << "}" << endl;

    ofs << padding << "}," << endl;
}

string Stats::list2string(double* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string Stats::list2string(double* list, int size, long* coords) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        // coords is 1,2,3,...
        long start = 0;
        if(i>0)
            start = coords[i-1];
        long end = coords[i];

        double total = 0.0;
        for(int k=start; k<end; k++)
            total += list[k];

        // get average
        if(end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string Stats::list2string(long* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

void Stats::reportHtml(ofstream& ofs, string filteringType, string readName) {
    reportHtmlQuality(ofs, filteringType, readName);
    reportHtmlContents(ofs, filteringType, readName);
}

bool Stats::isLongRead() {
    return mCycles > 300;
}

void Stats::reportHtmlQuality(ofstream& ofs, string filteringType, string readName) {

    // quality
    string subsection = filteringType + ": " + readName + ": quality";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    string alphabets[5] = {"A", "T", "C", "G", "mean"};
    string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<5; b++) {
        string base = alphabets[b];
        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mQualityCurves[base], total, x) + "],";
        json_str += "name: '" + base + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'quality'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

void Stats::reportHtmlContents(ofstream& ofs, string filteringType, string readName) {

    // content
    string subsection = filteringType + ": " + readName + ": base contents";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<6; b++) {
        string base = alphabets[b];
        long count = 0;
        if(base.size()==1) {
            char b = base[0] & 0x07;
            count = mBaseContents[b];
        } else {
            count = mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07] ;
        }
        string percentage = to_string((double)count * 100.0 / mBases);
        if(percentage.length()>5)
            percentage = percentage.substr(0,5);
        string name = base + "(" + percentage + "%)"; 

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mContentCurves[base], total, x) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'base content ratios'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

Stats* Stats::merge(vector<Stats*>& list) {
    if(list.size() == 0)
        return NULL;

    //get the most long cycles
    int cycles = 0;
    for(int t=0; t<list.size(); t++) {
        list[t]->summarize();
        cycles = max(cycles, list[t]->getCycles());
    }

    Stats* s = new Stats(list[0]->mOptions, list[0]->mIsRead2, cycles, 0);

    // init overrepresented seq maps
    map<string, long>::iterator iter;

    for(int t=0; t<list.size(); t++) {
        int curCycles =  list[t]->getCycles();
        // merge read number
        s->mReads += list[t]->mReads;
        s->mLengthSum += list[t]->mLengthSum;

        // merge per cycle counting for different bases
        for(int i=0; i<8; i++){
            for(int j=0; j<cycles && j<curCycles; j++) {
                s->mCycleQ30Bases[i][j] += list[t]->mCycleQ30Bases[i][j];
                s->mCycleQ20Bases[i][j] += list[t]->mCycleQ20Bases[i][j];
                s->mCycleBaseContents[i][j] += list[t]->mCycleBaseContents[i][j];
                s->mCycleBaseQual[i][j] += list[t]->mCycleBaseQual[i][j];
            }
        }

        // merge per cycle counting for all bases
        for(int j=0; j<cycles && j<curCycles; j++) {
            s->mCycleTotalBase[j] += list[t]->mCycleTotalBase[j];
            s->mCycleTotalQual[j] += list[t]->mCycleTotalQual[j];
        }

        /**最终mkmer没有使用到
        // merge kMer
        for(int i=0; i<s->mKmerBufLen; i++) {
            s->mKmer[i] += list[t]->mKmer[i];
        }
        **/
    }

    s->summarize();

    return s;
}

