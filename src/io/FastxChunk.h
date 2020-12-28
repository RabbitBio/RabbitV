/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00

  last modified by Zekun Yin 2020/5/18
*/

#ifndef H_FASTX_CHUNK
#define H_FASTX_CHUNK

#include <iostream>
#include <vector>

#include "Buffer.h"
#include "Common.h"
#include "DataPool.h"
#include "DataQueue.h"
#include "Globals.h"
#include "utils.h"

namespace mash {

namespace fa {

typedef core::DataChunk FastaDataChunk;

typedef core::TDataQueue<FastaDataChunk> FastaDataQueue;
typedef core::TDataPool<FastaDataChunk> FastaDataPool;

struct FastaChunk {
  FastaDataChunk *chunk;
  uint64 start;
  uint64 end;
  uint64 nseqs;

  void print() {
    std::cout << "chunk start: " << this->start << std::endl;
    std::cout << "chunk end: " << this->end << std::endl;
    std::cout << "chunk nseqs: " << this->nseqs << std::endl;
    return;
  }
};

}  // namespace fa

namespace fq {

typedef core::DataChunk FastqDataChunk;
typedef core::DataPairChunk FastqDataPairChunk;

typedef core::TDataQueue<FastqDataChunk> FastqDataQueue;
typedef core::TDataPool<FastqDataChunk> FastqDataPool;

struct FastqChunk {
  FastqDataChunk *chunk;
};

struct FastqPairChunk {
  FastqDataPairChunk *chunk;
};

}  // namespace fq

}  // namespace mash

#endif
