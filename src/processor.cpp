#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
}


Processor::~Processor(){
}

bool Processor::process() {
    //processing according whether the input data is pair-end or single-end.
    //for pair-end and single-end data, their quality control operations are different.
    if(mOptions->isPaired()) {
        PairEndProcessor p(mOptions);
        p.process();
    } else {
        SingleEndProcessor p(mOptions);
        p.process();
    }

    return true;
}