#ifndef __REFERENCE_H_
#define __REFERENCE_H_

#include <string>
#include <vector>
#include "Globals.h"
struct Reference {
  std::string name;
  std::string comment;
  std::string seq;
  std::string quality;
  std::string strand;
  uint64_t length;
  uint64_t gid;
};

struct neoReference {
  uint64_t pname;  // name offset
  uint64_t pcom;   // comment ??
  uint64_t pseq;   // sequence
  uint64_t pqual;  // quality
  uint64_t pstrand;

  uint64_t lname;
  uint64_t lcom;
  uint64_t lseq;
  uint64_t lqual;
  uint64_t lstrand;
  // uint64_t length;
  uint64_t gid;
  mash::byte *base;
};
typedef std::vector<Reference> SeqInfos;
typedef Reference OneSeqInfo;

#endif