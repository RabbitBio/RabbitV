#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

using namespace std;

#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

class DuplicationOptions {
public:
    DuplicationOptions() {
        enabled = true;
        keylen = 12;
        histSize = 32;
    }
public:
    bool enabled;
    int keylen;
    int histSize;
};


class LowComplexityFilterOptions {
public:
    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
public:
    bool enabled;
    double threshold;
};

class PolyGTrimmerOptions {
public:
    PolyGTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class PolyXTrimmerOptions {
public:
    PolyXTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class UMIOptions {
public:
    UMIOptions() {
        enabled = false;
        location = UMI_LOC_NONE;
        length = 0;
        skip = 0;
    }
public:
    bool enabled;
    int location;
    int length;
    int skip;
    string prefix;
    string separator;
};

class CorrectionOptions {
public:
    CorrectionOptions() {
        enabled = false;
    }
public:
    bool enabled;
};

class QualityCutOptions {
public:
    QualityCutOptions() {
        enabledFront = false;
        enabledTail = false;
        enabledRight = false;
        windowSizeShared = 4;
        qualityShared = 20;
        windowSizeFront = windowSizeShared;
        qualityFront = qualityShared;
        windowSizeTail = windowSizeShared;
        qualityTail = qualityShared;
        windowSizeRight = windowSizeShared;
        qualityRight = qualityShared;
    }
public:
    // enable 5' cutting by quality
    bool enabledFront;
    // enable 3' cutting by quality
    bool enabledTail;
    // enable agressive cutting mode
    bool enabledRight;
    // the sliding window size
    int windowSizeShared;
    // the mean quality requirement
    int qualityShared;
    // the sliding window size for cutting by quality in 5'
    int windowSizeFront;
    // the mean quality requirement for cutting by quality in 5'
    int qualityFront;
    // the sliding window size for cutting by quality in 3'
    int windowSizeTail;
    // the mean quality requirement for cutting by quality in 3'
    int qualityTail;
    // the sliding window size for cutting by quality in aggressive mode
    int windowSizeRight;
    // the mean quality requirement for cutting by quality in aggressive mode
    int qualityRight;
};

class AdapterOptions {
public:
    AdapterOptions() {
        enabled = true;
        hasSeqR1 = false;
        hasSeqR2 = false;
        detectAdapterForPE = false;
    }
public:
    bool enabled;
    string sequence;
    string sequenceR2;
    string detectedAdapter1;
    string detectedAdapter2;
    vector<string> seqsInFasta;
    string fastaFile;
    bool hasSeqR1;
    bool hasSeqR2;
    bool hasFasta;
    bool detectAdapterForPE;
};

class TrimmingOptions {
public:
    TrimmingOptions() {
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
        maxLen1 = 0;
        maxLen2 = 0;
    }
public:
    // trimming first cycles for read1
    int front1;
    // trimming last cycles for read1
    int tail1;
    // trimming first cycles for read2
    int front2;
    // trimming last cycles for read2
    int tail2;
    // max length of read1
    int maxLen1;
    // max length of read2
    int maxLen2;
};

class QualityFilteringOptions {
public:
    QualityFilteringOptions() {
        enabled = true;
        // '0' = Q15
        qualifiedQual = '0';
        unqualifiedPercentLimit = 40;
        nBaseLimit = 5;
    }
public:
    // quality filter enabled
    bool enabled;
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    char qualifiedQual;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedPercentLimit;
    // if n_base_number > nBaseLimit, then discard this read
    int nBaseLimit;
    // if average qual score < avgQualReq, then discard this read
    int avgQualReq;
};

class ReadLengthFilteringOptions {
public:
    ReadLengthFilteringOptions() {
        enabled = false;
        requiredLength = 15;
        maxLength = 0;
    }
public:
    // length filter enabled
    bool enabled;
    // if read_length < requiredLength, then this read is discard
    int requiredLength;
    // length limit, 0 for no limitation
    int maxLength;
};

class Options{
public:
    Options();
    void init();
    bool isPaired();
    bool validate();
    bool adapterCuttingEnabled();
    bool polyXTrimmingEnabled();
    string getAdapter1();
    string getAdapter2();
    bool shallDetectAdapter(bool isR2 = false);
    void loadFastaAdapters();

public:
    // file name of read1 input
    string in1;
    // file name of read2 input
    string in2;
    // file name of read1 output
    string out1;
    // file name of read2 output
    string out2;
    // genome FASTA file
    string genomeFile;
    // kmer FASTA file
    string kmerFile;
    // kmer FASTA file
    string kmerCollectionFile;
    // json file
    string jsonFile;
    // html file
    string htmlFile;
    // html report title
    string reportTitle;
    // compression level
    int compression;
    // the input file is using phred64 quality scoring
    bool phred64;
    // do not rewrite existing files
    bool dontOverwrite;
    // read STDIN
    bool inputFromSTDIN;
    // write STDOUT
    bool outputToSTDOUT;
    // the input R1 file is interleaved
    bool interleavedInput;
    // only process first N reads
    int readsToProcess;
    // worker thread number
    int thread;
    // trimming options
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;
    // length filtering options
    ReadLengthFilteringOptions lengthFilter;
    // adapter options
    AdapterOptions adapter;
    // options for quality cutting
    QualityCutOptions qualityCut;
    // options for base correction
    CorrectionOptions correction;
    // options for UMI
    UMIOptions umi;
    // 3' end polyG trimming, default for Illumina NextSeq/NovaSeq
    PolyGTrimmerOptions polyGTrim;
    // 3' end polyX trimming
    PolyXTrimmerOptions polyXTrim;
    int seqLen1;
    int seqLen2;
    // low complexity filtering
    LowComplexityFilterOptions complexityFilter;
    // options for duplication profiling
    DuplicationOptions duplicate;
    // options for duplication profiling
    int insertSizeMax;
    // overlap analysis threshold
    int overlapRequire;
    int overlapDiffLimit;
    int overlapDiffPercentLimit;
    // output debug information
    bool verbose;
    // the length of KMER, default is 25
    int kmerKeyLen;
    // the threshold of positive result
    double positiveThreshold;
    // the threshold of depth for a region considered as covered
    double depthThreshold;
    // if ed(read, genome) <= edThreshold, then think it as a match
    int edThreshold;
    // the bin size for stat coverage and edit distance
    int statsBinSize;
    // read with length >= longReadThreshold will be considered as long reads
    int longReadThreshold;
    // long reads will be split to reads with length <= segmentLength
    int segmentLength;
    //  coverage threshold to be reported in kmer collection results
    double kcCoverageThreshold;
    // coverage for high-confidence KCR
    double kcCoverageHighConfidence;
    // median hit for high-confidence KCR
    double kcMedianHitHighConfidence;

};

#endif