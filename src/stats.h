#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "read.h"
#include "options.h"

using namespace std;

class Stats{
public:
    // this @guessedCycles parameter should be calculated using the first several records
    Stats(Options* opt, bool isRead2 = false, int guessedCycles = 0, int bufferMargin = 1024);
    ~Stats();
    int getCycles();
    long getReads();
    long getBases();
    long getQ20();
    long getQ30();
    long getGCNumber();
    // by default the qualified qual score is Q20 ('5')
    void statRead(Read* r);

    static Stats* merge(vector<Stats*>& list);
    void print();
    void summarize(bool forced = false);
    // a port of JSON report
    void reportJson(ofstream& ofs, string padding);
    // a port of HTML report
    void reportHtml(ofstream& ofs, string filteringType, string readName);
    void reportHtmlQuality(ofstream& ofs, string filteringType, string readName);
    void reportHtmlContents(ofstream& ofs, string filteringType, string readName);
    bool isLongRead();
    int getMeanLength();

public:
    static string list2string(double* list, int size);
    static string list2string(double* list, int size, long* coords);
    static string list2string(long* list, int size);
    static int base2val(char base);

private:
    void extendBuffer(int newBufLen);
    void countmCycleTotalBase();//利用mCycleBaseContents来总和得到mCycleTotalBase
    void countmCycleTotalQual();//利用mCycleBaseQual来总和得到mCycleTotalQual
    void moveTempArr2FinalArr();

private:
    Options* mOptions;
    bool mIsRead2;
    long mReads;
    int mEvaluatedSeqLen;
    /* 
    why we use 8 here?
    map A/T/C/G/N to 0~7 by their ASCII % 8:
    'A' % 8 = 1
    'T' % 8 = 4
    'C' % 8 = 3
    'G' % 8 = 7
    'N' % 8 = 6
    */
    long *mCycleQ30Bases[8];
    long *mCycleQ20Bases[8];
    long *mCycleBaseContents[8];
    long *mCycleBaseQual[8];
    long *mCycleTotalBase;
    long *mCycleTotalQual;

    /**
     * 因为使用了向量化，由于下标是用int表示的，所以如果使用int的话，会比使用long每次要多进行一倍的运算
     * 所以采用int大小的临时数组
     **/
    int *mCycleQ30BasesTemp;
    int *mCycleQ20BasesTemp;
    int *mCycleBaseContentsTemp;
    int *mCycleBaseQualTemp;

    //记录临时数组还能容纳多少条序列
    int tempEmpty;

    //一个int最多记录2^31-1个，由于fastq序列的碱基质量范围在0-40，所以最多(2^31-1)/40条序列就可能会满
    //这个时候就需要从临时数组转移到实际数组了，而(2^31-1)/40=53687091，这里取50000000满足要求
    static const int maxCount = 50000000;
    /**因为这玩意并没有对程序其他部分还有最终需要的输出结果产生影响，所以可以直接删去
    long *mKmer;
    **/
   
    map<string, double*> mQualityCurves;
    map<string, double*> mContentCurves;


    int mCycles;
    int mBufLen;
    long mBases;
    long mQ20Bases[8];
    long mQ30Bases[8];
    long mBaseContents[8];
    long mQ20Total;
    long mQ30Total;
    bool summarized;
    long mKmerMax;
    long mKmerMin;
    int mKmerBufLen;
    long mLengthSum;
};

#endif