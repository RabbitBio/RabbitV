#include "jsonreporter.h"

JsonReporter::JsonReporter(Options* opt){
    mOptions = opt;
    mDupHist = NULL;
    mDupRate = 0;
}

JsonReporter::~JsonReporter(){
}

void JsonReporter::setDupHist(int* dupHist, double* dupMeanGC, double dupRate) {
    mDupHist = dupHist;
    mDupMeanGC = dupMeanGC;
    mDupRate = dupRate;
}

void JsonReporter::setInsertHist(long* insertHist, int insertSizePeak) {
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

extern string command;
void JsonReporter::report(VirusDetector* vd, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    ofstream ofs;
    ofs.open(mOptions->jsonFile, ifstream::out);
    ofs << "{" << endl;

    long pre_total_reads = preStats1->getReads();
    if(preStats2)
        pre_total_reads += preStats2->getReads();

    long pre_total_bases = preStats1->getBases();
    if(preStats2)
        pre_total_bases += preStats2->getBases();

    long pre_q20_bases = preStats1->getQ20();
    if(preStats2)
        pre_q20_bases += preStats2->getQ20();

    long pre_q30_bases = preStats1->getQ30();
    if(preStats2)
        pre_q30_bases += preStats2->getQ30();

    long pre_total_gc = preStats1->getGCNumber();
    if(preStats2)
        pre_total_gc += preStats2->getGCNumber();

    long post_total_reads = postStats1->getReads();
    if(postStats2)
        post_total_reads += postStats2->getReads();

    long post_total_bases = postStats1->getBases();
    if(postStats2)
        post_total_bases += postStats2->getBases();

    long post_q20_bases = postStats1->getQ20();
    if(postStats2)
        post_q20_bases += postStats2->getQ20();

    long post_q30_bases = postStats1->getQ30();
    if(postStats2)
        post_q30_bases += postStats2->getQ30();

    long post_total_gc = postStats1->getGCNumber();
    if(postStats2)
        post_total_gc += postStats2->getGCNumber();

    // KMER detection
    Kmer* kmer = vd->getKmer();
    if(kmer) {
        string detectionResult;
        if(kmer->getMeanHit() >= mOptions->positiveThreshold)
            detectionResult = "POSITIVE";
        else
            detectionResult = "NEGATIVE";

        ofs << "\t" << "\"kmer_detection_result\": {" << endl;
        ofs << "\t\t" << "\"result\": \"" << detectionResult << "\"," << endl;
        ofs << "\t\t" << "\"mean_coverage\": " << kmer->getMeanHit() << "," << endl;
        ofs << "\t\t" << "\"positive_thread\": " << mOptions->positiveThreshold << "," << endl;

        // unique kmer hits
        ofs << "\t\t" << "\"kmer_hits\": {" << endl;
            kmer->reportJSON(ofs);
        ofs << "\t\t" << "}" << endl;


        ofs << "\t" << "}," << endl;
    }

    // KMER detection
    Genomes* genome = vd->getGenomes();
    if(genome) {
        genome->reportJSON(ofs);
    }

    // KMER detection
    KmerCollection* kc = vd->getKmerCollection();
    if(kc) {
        kc->reportJSON(ofs);
    }

    // summary
    ofs << "\t" << "\"summary\": {" << endl;

    ofs << "\t\t" << "\"before_filtering\": {" << endl;
    ofs << "\t\t\t" << "\"total_reads\":" << pre_total_reads << "," << endl; 
    ofs << "\t\t\t" << "\"total_bases\":" << pre_total_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_bases\":" << pre_q20_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q30_bases\":" << pre_q30_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_rate\":" << (pre_total_bases == 0?0.0:(double)pre_q20_bases / (double)pre_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"q30_rate\":" << (pre_total_bases == 0?0.0:(double)pre_q30_bases / (double)pre_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"read1_mean_length\":" << preStats1->getMeanLength() << "," << endl;
    if(mOptions->isPaired())
        ofs << "\t\t\t" << "\"read2_mean_length\":" << preStats2->getMeanLength() << "," << endl;
    ofs << "\t\t\t" << "\"gc_content\":" << (pre_total_bases == 0?0.0:(double)pre_total_gc / (double)pre_total_bases)  << endl; 
    ofs << "\t\t" << "}," << endl;

    ofs << "\t\t" << "\"after_filtering\": {" << endl;
    ofs << "\t\t\t" << "\"total_reads\":" << post_total_reads << "," << endl; 
    ofs << "\t\t\t" << "\"total_bases\":" << post_total_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_bases\":" << post_q20_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q30_bases\":" << post_q30_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_rate\":" << (post_total_bases == 0?0.0:(double)post_q20_bases / (double)post_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"q30_rate\":" << (post_total_bases == 0?0.0:(double)post_q30_bases / (double)post_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"read1_mean_length\":" << postStats1->getMeanLength() << "," << endl;
    if(mOptions->isPaired())
        ofs << "\t\t\t" << "\"read2_mean_length\":" << postStats2->getMeanLength() << "," << endl;
    ofs << "\t\t\t" << "\"gc_content\":" << (post_total_bases == 0?0.0:(double)post_total_gc / (double)post_total_bases)  << endl; 
    ofs << "\t\t" << "}";

    ofs << endl;

    ofs << "\t" << "}," << endl;

    if(result) {
        ofs << "\t" << "\"filtering_result\": " ;
        result -> reportJson(ofs, "\t");
    }

    if(mOptions->duplicate.enabled) {
        ofs << "\t" << "\"duplication\": {" << endl;
        ofs << "\t\t\"rate\": " << mDupRate << "," << endl;
        ofs << "\t\t\"histogram\": [";
        for(int d=1; d<mOptions->duplicate.histSize; d++) {
            ofs << mDupHist[d];
            if(d!=mOptions->duplicate.histSize-1)
                ofs << ",";
        }
        ofs << "]," << endl;
        ofs << "\t\t\"mean_gc\": [";
        for(int d=1; d<mOptions->duplicate.histSize; d++) {
            ofs << mDupMeanGC[d];
            if(d!=mOptions->duplicate.histSize-1)
                ofs << ",";
        }
        ofs << "]" << endl;
        ofs << "\t" << "}";
        ofs << "," << endl;
    }

    if(mOptions->isPaired()) {
        ofs << "\t" << "\"insert_size\": {" << endl;
        ofs << "\t\t\"peak\": " << mInsertSizePeak << "," << endl;
        ofs << "\t\t\"unknown\": " << mInsertHist[mOptions->insertSizeMax] << "," << endl;
        ofs << "\t\t\"histogram\": [";
        for(int d=0; d<mOptions->insertSizeMax; d++) {
            ofs << mInsertHist[d];
            if(d!=mOptions->insertSizeMax-1)
                ofs << ",";
        }
        ofs << "]" << endl;
        ofs << "\t" << "}";
        ofs << "," << endl;
    }

    if(result && mOptions->adapterCuttingEnabled()) {
        ofs << "\t" << "\"adapter_cutting\": " ;
        result -> reportAdapterJson(ofs, "\t");
    }

    if(result && mOptions->polyXTrimmingEnabled()) {
        ofs << "\t" << "\"polyx_trimming\": " ;
        result -> reportPolyXTrimJson(ofs, "\t");
    }

    if(preStats1) {
        ofs << "\t" << "\"read1_before_filtering\": " ;
        preStats1 -> reportJson(ofs, "\t");
    }

    if(preStats2) {
        ofs << "\t" << "\"read2_before_filtering\": " ;
        preStats2 -> reportJson(ofs, "\t");
    }

    if(postStats1) {
        string name = "read1_after_filtering";
        ofs << "\t" << "\"" << name << "\": " ;
        postStats1 -> reportJson(ofs, "\t");
    }

    if(postStats2) {
        ofs << "\t" << "\"read2_after_filtering\": " ;
        postStats2 -> reportJson(ofs, "\t");
    }

    ofs << "\t\"command\": " << "\"" << command << "\"" << endl;

    ofs << "}";
}