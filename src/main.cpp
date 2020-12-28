#include <stdio.h>
#include "fastqreader.h"
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"

// TODO: code refactoring to remove these global variables
string command;
mutex logmtx;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cerr << "fastv: an ultra-fast tool for fast identification of SARS-CoV-2 and other microbes from sequencing data." << endl << "version " << FASTV_VER << endl;
        //cerr << "fastv --help to see the help"<<endl;
        //return 0;
    }
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "fastv " << FASTV_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    cmd.add<string>("in1", 'i', "read1 input file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("out1", 'o', "file name to store read1 with on-target sequences", false, "");
    cmd.add<string>("out2", 'O', "file name to store read2 with on-target sequences", false, "");
    cmd.add<string>("kmer_collection", 'c', "the unique k-mer collection file in fasta format, see an example: http://opengene.org/kmer_collection.fasta", false, "");
    cmd.add<string>("kmer", 'k', "the unique k-mer file of the detection target in fasta format. data/SARS-CoV-2.kmer.fa will be used if none of k-mer/Genomes/k-mer_Collection file is specified", false, "");
    cmd.add<string>("genomes", 'g', "the genomes file of the detection target in fasta format. data/SARS-CoV-2.genomes.fa will be used if none of k-mer/Genomes/k-mer_Collection file is specified", false, "");
    cmd.add<float>("positive_threshold", 'p', "the data is considered as POSITIVE, when its mean coverage of unique kmer >= positive_threshold (0.001 ~ 100). 0.1 by default.", false, 0.1);
    cmd.add<float>("depth_threshold", 'd', "For coverage calculation. A region is considered covered when its mean depth >= depth_threshold (0.001 ~ 1000). 1.0 by default.", false, 1.0);
    cmd.add<int>("ed_threshold", 'E', "If the edit distance of a sequence and a genome region is <=ed_threshold, then consider it a match (0 ~ 50). 8 by default.", false, 8);
    cmd.add<int>("long_read_threshold", 0, "A read will be considered as long read if its length >= long_read_threshold (100 ~ 10000). 200 by default.", false, 200);
    cmd.add<int>("read_segment_len", 0, "A long read will be splitted to read segments, with each <= read_segment_len (50 ~ 5000, should be < long_read_threshold). 100 by default.", false, 100);
    cmd.add<int>("bin_size", 0, "For coverage calculation. The genome is splitted to many bins, with each bin has a length of bin_size (1 ~ 100000), default 0 means adaptive.", false, 0);
    cmd.add<double>("kc_coverage_threshold", 0, "For each genome in the k-mer collection FASTA, report it when its coverage > kc_coverage_threshold. Default is 0.01.", false, 0.01);
    cmd.add<double>("kc_high_confidence_coverage_threshold", 0, "For each genome in the k-mer collection FASTA, report it as high confidence when its coverage > kc_high_confidence_coverage_threshold. Default is 0.9.", false, 0.9);
    cmd.add<int>("kc_high_confidence_median_hit_threshold", 0, "For each genome in the k-mer collection FASTA, report it as high confidence when its median hits > kc_high_confidence_median_hit_threshold. Default is 5.", false, 5);

    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "fastv.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "fastv.html");
    cmd.add<string>("report_title", 'R', "should be quoted with \' or \", default is \"fastv report\"", false, "fastv report");

    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 4", false, 4);

    // qother I/O
    cmd.add("phred64", '6', "indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
    cmd.add("stdin", 0, "input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.");
    cmd.add("stdout", 0, "stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.");
    cmd.add("interleaved_in", 0, "indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.");
    cmd.add<int>("reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.", false, 0);
    cmd.add("dont_overwrite", 0, "don't overwrite existing files. Overwritting is allowed by default.");
    cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");

    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.", false, "auto");
    cmd.add<string>("adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>", false, "auto");
    cmd.add<string>("adapter_fasta", 0, "specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file", false, "");
    cmd.add("detect_adapter_for_pe", 0, "by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.");

    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("max_len1", 'b', "if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("max_len2", 'B', "if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings", false, 0);

    // polyG tail trimming
    cmd.add<int>("poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default.", false, 10);
    cmd.add("disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    
    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);

    // cutting by quality
    cmd.add("cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_right", 'r', "move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.");
    cmd.add<int>("cut_window_size", 'W', "the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)", false, 20);
    cmd.add<int>("cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_right_window_size", 0, "the window size option of cut_right, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_right_mean_quality", 0, "the mean quality requirement option for cut_right, default to cut_mean_quality if not specified", false, 20);


    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.", false, 15);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);
    cmd.add<int>("average_qual", 'e', "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement", false, 0);

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 15.", false, 15);
    cmd.add<int>("length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.", false, 0);

    // low complexity filtering
    cmd.add("low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);

    // filter by indexes
    cmd.add<string>("filter_by_index1", 0, "specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line", false, "");
    cmd.add<string>("filter_by_index2", 0, "specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line", false, "");
    cmd.add<int>("filter_by_index_threshold", 0, "the allowed difference of index barcode for index filtering, default 0 means completely identical.", false, 0);
    
    // base correction in overlapped regions of paired end data
    cmd.add("correction", 'C', "enable base correction in overlapped regions (only for PE data), default is disabled");
    cmd.add<int>("overlap_len_require", 0, "the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default.", false, 30);
    cmd.add<int>("overlap_diff_limit", 0, "the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.", false, 5);
    cmd.add<int>("overlap_diff_percent_limit", 0, "the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.", false, 20);

    // umi
    cmd.add("umi", 'U', "enable unique molecular identifier (UMI) preprocessing");
    cmd.add<string>("umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none", false, "");
    cmd.add<int>("umi_len", 0, "if the UMI is in read1/read2, its length should be provided", false, 0);
    cmd.add<string>("umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default", false, "");
    cmd.add<int>("umi_skip", 0, "if the UMI is in read1/read2, fastv can skip several bases following UMI, default is 0", false, 0);

    
    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
        return 0;
    }

    Options opt;

    // I/O
    opt.in1 = cmd.get<string>("in1");
    opt.in2 = cmd.get<string>("in2");
    opt.out1 = cmd.get<string>("out1");
    opt.out2 = cmd.get<string>("out2");

    string fastvProgPath = string(argv[0]);
    string fastvDir = dirname(fastvProgPath);
    opt.kmerFile = cmd.get<string>("kmer");
    opt.kmerCollectionFile = cmd.get<string>("kmer_collection");
    opt.genomeFile = cmd.get<string>("genomes");
    if(opt.kmerFile.empty() && opt.genomeFile.empty() && opt.kmerCollectionFile.empty()) {
        cerr << endl << "SARS-CoV-2 Detection Mode..." << endl;
        cerr << "Since none of k-mer file (-k), Genomes file (-g) and k-mer_Collection file (-c) is specified, fastv will try to load SARS-CoV-2 k-mer/Genomes files from " << joinpath(fastvDir, "data") << endl;
        string kmerFile = joinpath(fastvDir, "data/SARS-CoV-2.kmer.fa");
        if(file_exists(kmerFile)) {
            cerr << "Found k-mer file: " << kmerFile << endl;
            opt.kmerFile = kmerFile;
        } else {
            cerr << "Didn't find k-mer file: " << kmerFile << endl;
        }
        string genomeFile = joinpath(fastvDir, "data/SARS-CoV-2.genomes.fa");
        if(file_exists(genomeFile)) {
            cerr << "Found Genomes file: " << genomeFile << endl;
            opt.genomeFile = genomeFile;
        } else {
            cerr << "Didn't find Genomes file: " << genomeFile << endl;
        }

        if(!file_exists(kmerFile) && !file_exists(genomeFile)) {
            cerr << "Could't find the built-in k-mer file or Genomes file " << endl;
            error_exit("Please specify at least one k-mer file (-k) or one Genomes file (-g)"); 
        }
        cerr << endl;
    }

    opt.positiveThreshold = cmd.get<float>("positive_threshold");
    opt.depthThreshold = cmd.get<float>("depth_threshold");
    opt.edThreshold = cmd.get<int>("ed_threshold");
    opt.longReadThreshold = cmd.get<int>("long_read_threshold");
    opt.segmentLength = cmd.get<int>("read_segment_len");
    opt.statsBinSize = cmd.get<int>("bin_size");

    opt.kcCoverageThreshold = cmd.get<double>("kc_coverage_threshold");
    opt.kcCoverageHighConfidence = cmd.get<double>("kc_high_confidence_coverage_threshold");
    opt.kcMedianHitHighConfidence = cmd.get<int>("kc_high_confidence_median_hit_threshold");

    opt.compression = cmd.get<int>("compression");
    opt.readsToProcess = cmd.get<int>("reads_to_process");
    opt.phred64 = cmd.exist("phred64");
    opt.dontOverwrite = cmd.exist("dont_overwrite");
    opt.inputFromSTDIN = cmd.exist("stdin");
    opt.outputToSTDOUT = cmd.exist("stdout");
    opt.interleavedInput = cmd.exist("interleaved_in");
    opt.verbose = cmd.exist("verbose");

    // adapter cutting
    opt.adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt.adapter.detectAdapterForPE = cmd.exist("detect_adapter_for_pe");
    opt.adapter.sequence = cmd.get<string>("adapter_sequence");
    opt.adapter.sequenceR2 = cmd.get<string>("adapter_sequence_r2");
    opt.adapter.fastaFile = cmd.get<string>("adapter_fasta");
    if(opt.adapter.sequenceR2=="auto" && !opt.adapter.detectAdapterForPE && opt.adapter.sequence != "auto") {
        opt.adapter.sequenceR2 = opt.adapter.sequence;
    }
    if(!opt.adapter.fastaFile.empty()) {
        opt.loadFastaAdapters();
    }

    // trimming
    opt.trim.front1 = cmd.get<int>("trim_front1");
    opt.trim.tail1 = cmd.get<int>("trim_tail1");
    opt.trim.maxLen1 = cmd.get<int>("max_len1");
    // read2 settings follows read1 if it's not specified
    if(cmd.exist("trim_front2"))
        opt.trim.front2 = cmd.get<int>("trim_front2");
    else
        opt.trim.front2 = opt.trim.front1;
    if(cmd.exist("trim_tail2"))
        opt.trim.tail2 = cmd.get<int>("trim_tail2");
    else
        opt.trim.tail2 = opt.trim.tail1;
    if(cmd.exist("max_len2"))
        opt.trim.maxLen2 = cmd.get<int>("max_len2");
    else
        opt.trim.maxLen2 = opt.trim.maxLen1;

    // polyG tail trimming
    if(cmd.exist("disable_trim_poly_g")) {
        opt.polyGTrim.enabled = false;
    }
    opt.polyGTrim.minLen = cmd.get<int>("poly_g_min_len");

    // polyX tail trimming
    if(cmd.exist("trim_poly_x")) {
        opt.polyXTrim.enabled = true;
    }
    opt.polyXTrim.minLen = cmd.get<int>("poly_x_min_len");


    // sliding window cutting by quality
    opt.qualityCut.enabledFront = cmd.exist("cut_front");
    opt.qualityCut.enabledTail = cmd.exist("cut_tail");
    opt.qualityCut.enabledRight = cmd.exist("cut_right");

    opt.qualityCut.windowSizeShared = cmd.get<int>("cut_window_size");
    opt.qualityCut.qualityShared = cmd.get<int>("cut_mean_quality");

    if(cmd.exist("cut_front_window_size"))
        opt.qualityCut.windowSizeFront = cmd.get<int>("cut_front_window_size");
    else
        opt.qualityCut.windowSizeFront = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_front_mean_quality"))
        opt.qualityCut.qualityFront = cmd.get<int>("cut_front_mean_quality");
    else
        opt.qualityCut.qualityFront = opt.qualityCut.qualityShared;

    if(cmd.exist("cut_tail_window_size"))
        opt.qualityCut.windowSizeTail = cmd.get<int>("cut_tail_window_size");
    else
        opt.qualityCut.windowSizeTail = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_tail_mean_quality"))
        opt.qualityCut.qualityTail = cmd.get<int>("cut_tail_mean_quality");
    else
        opt.qualityCut.qualityTail = opt.qualityCut.qualityShared;

    if(cmd.exist("cut_right_window_size"))
        opt.qualityCut.windowSizeRight = cmd.get<int>("cut_right_window_size");
    else
        opt.qualityCut.windowSizeRight = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_right_mean_quality"))
        opt.qualityCut.qualityRight = cmd.get<int>("cut_right_mean_quality");
    else
        opt.qualityCut.qualityRight = opt.qualityCut.qualityShared;

    // raise a warning if cutting option is not enabled but -W/-M is enabled
    if(!opt.qualityCut.enabledFront && !opt.qualityCut.enabledTail && !opt.qualityCut.enabledRight) {
        if(cmd.exist("cut_window_size") || cmd.exist("cut_mean_quality") 
            || cmd.exist("cut_front_window_size") || cmd.exist("cut_front_mean_quality") 
            || cmd.exist("cut_tail_window_size") || cmd.exist("cut_tail_mean_quality") 
            || cmd.exist("cut_right_window_size") || cmd.exist("cut_right_mean_quality"))
            cerr << "WARNING: you specified the options for cutting by quality, but forogt to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    }

    // quality filtering
    opt.qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt.qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt.qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt.qualfilter.avgQualReq = cmd.get<int>("average_qual");
    opt.qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt.lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt.lengthFilter.requiredLength = cmd.get<int>("length_required");
    opt.lengthFilter.maxLength = cmd.get<int>("length_limit");

    // low complexity filter
    opt.complexityFilter.enabled = cmd.exist("low_complexity_filter");
    opt.complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;

    // overlap correction
    opt.correction.enabled = cmd.exist("correction");
    opt.overlapRequire = cmd.get<int>("overlap_len_require");
    opt.overlapDiffLimit = cmd.get<int>("overlap_diff_limit");
    opt.overlapDiffPercentLimit = cmd.get<int>("overlap_diff_percent_limit");

    // threading
    opt.thread = cmd.get<int>("thread");

    // reporting
    opt.jsonFile = cmd.get<string>("json");
    opt.htmlFile = cmd.get<string>("html");
    opt.reportTitle = cmd.get<string>("report_title");

    // umi
    opt.umi.enabled = cmd.exist("umi");
    opt.umi.length = cmd.get<int>("umi_len");
    opt.umi.prefix = cmd.get<string>("umi_prefix");
    opt.umi.skip = cmd.get<int>("umi_skip");
    if(opt.umi.enabled) {
        string umiLoc = cmd.get<string>("umi_loc");
        str2lower(umiLoc);
        if(umiLoc.empty())
            error_exit("You've enabled UMI by (--umi), you should specify the UMI location by (--umi_loc)");
        if(umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" && umiLoc != "per_index" && umiLoc != "per_read") {
            error_exit("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }
        if(!opt.isPaired() && (umiLoc == "index2" || umiLoc == "read2"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if(opt.umi.length == 0 && (umiLoc == "read1" || umiLoc == "read2" ||  umiLoc == "per_read"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the length is not specified (--umi_len).");
        if(umiLoc == "index1") {
            opt.umi.location = UMI_LOC_INDEX1;
        } else if(umiLoc == "index2") {
            opt.umi.location = UMI_LOC_INDEX2;
        } else if(umiLoc == "read1") {
            opt.umi.location = UMI_LOC_READ1;
        } else if(umiLoc == "read2") {
            opt.umi.location = UMI_LOC_READ2;
        } else if(umiLoc == "per_index") {
            opt.umi.location = UMI_LOC_PER_INDEX;
        } else if(umiLoc == "per_read") {
            opt.umi.location = UMI_LOC_PER_READ;
        }
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    bool supportEvaluation = !opt.inputFromSTDIN && opt.in1!="/dev/stdin";

    Evaluator eva(&opt);
    if(supportEvaluation) {
        eva.evaluateSeqLen();
    }

    long readNum = 0;

    // using evaluator to guess how many reads in total
    if(opt.shallDetectAdapter(false)) {
        if(!supportEvaluation) {
            //cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        }
        else {
            //cerr << "Detecting adapter sequence for read1..." << endl;
            string adapt = eva.evalAdapterAndReadNum(readNum, false);
            if(adapt.length() > 60 )
                adapt.resize(0, 60);
            if(adapt.length() > 0 ) {
                opt.adapter.sequence = adapt;
                opt.adapter.detectedAdapter1 = adapt;
            } else {
                //cerr << "No adapter detected for read1" << endl;
                opt.adapter.sequence = "";
            }
            cerr << endl;
        }
    }
    if(opt.shallDetectAdapter(true)) {
        if(!supportEvaluation) {
            //cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        }
        else {
            //cerr << "Detecting adapter sequence for read2..." << endl;
            string adapt = eva.evalAdapterAndReadNum(readNum, true);
            if(adapt.length() > 60 )
                adapt.resize(0, 60);
            if(adapt.length() > 0 ) {
                opt.adapter.sequenceR2 = adapt;
                opt.adapter.detectedAdapter2 = adapt;
            } else {
                //cerr << "No adapter detected for read2" << endl;
                opt.adapter.sequenceR2 = "";
            }
            cerr << endl;
        }
    }

    opt.validate();

    // using evaluator to check if it's two color system
    if(!cmd.exist("disable_trim_poly_g") && supportEvaluation) {
        bool twoColorSystem = eva.isTwoColorSystem();
        if(twoColorSystem){
            opt.polyGTrim.enabled = true;
        }
    }

    Processor p(&opt);
    p.process();
    
    time_t t2 = time(NULL);

    cerr << endl << "JSON report: " << opt.jsonFile << endl;
    cerr << "HTML report: " << opt.htmlFile << endl;
    cerr << endl << command << endl;
    cerr << "fastv v" << FASTV_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
