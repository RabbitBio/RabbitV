/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00

  last modified by Zekun Yin 2020/5/18
*/

#include <iostream>
#include <string>

#include "Reference.h"
#include "FastxStream.h"

namespace mash {

namespace fa {

FastaChunk *FastaFileReader::readNextChunk() {
  FastaDataChunk *part = NULL;
  recordsPool.Acquire(part);
  FastaChunk *dataPart = new FastaChunk;
  dataPart->chunk = part;
  if (ReadNextChunk(dataPart, this->seqInfos)) {
    return dataPart;
  } else {
    std::cout << "read over!" << std::endl;
    recordsPool.Release(part);
    return NULL;
  }
}

/*
Description:
  this function make sure one FastaChunk(dataPart) contains at least
  one whole sequence, because
*/
FastaChunk *FastaFileReader::readNextChunkList() {
  FastaDataChunk *part = NULL;
  recordsPool.Acquire(part);
  FastaChunk *dataPart = new FastaChunk;
  dataPart->chunk = part;
  FastaDataChunk *currnet = NULL;
  bool continue_read = false;
  if (ReadNextFaChunk(dataPart->chunk, this->seqInfos, continue_read)) {
    FastaDataChunk *current = part;
    while (continue_read) {
      FastaDataChunk *append = NULL;
      recordsPool.Acquire(append);
      if (ReadNextFaChunk(append, this->seqInfos, continue_read)) {
        current->next = append;
        currnet = append;
      } else {
        recordsPool.Release(append);
        break;
      }
    }
    return dataPart;
  } else {
    recordsPool.Release(part);
    return NULL;
  }
}

bool FastaFileReader::ReadNextChunk(FastaChunk *dataChunk_, SeqInfos &seqInfos) {
  FastaDataChunk *chunk_ = dataChunk_->chunk;
  if (Eof()) {
    chunk_->size = 0;
    return false;
  }

  // flush the data from previous incomplete chunk
  uchar *data = chunk_->data.Pointer();
  const uint64 cbufSize = chunk_->data.Size();
  chunk_->size = 0;
  int64 toRead = cbufSize - bufferSize;  // buffersize: size left from last chunk

  if (bufferSize > 0) {
    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
    chunk_->size = bufferSize;
    bufferSize = 0;
  }

  // read the next chunk
  int64 r = this->Read(data + chunk_->size, toRead);
  // std::cout << "r is :" << r << std::endl;
  // std::cout << "toRead: " << toRead << std::endl;

  if (r > 0) {
    if (r == toRead)  // somewhere before end
    {
      // dealing with halo region
      uint64 chunkEnd = FindCutPos(dataChunk_, data, cbufSize, mHalo, seqInfos);
      chunk_->size = chunkEnd;  // - 1; //1 char back from last '>'

      if (usesCrlf) chunk_->size -= 1;

      // copy tail to swapBuffer
      int haloCount = 0;
      uint64 tmpPtr = chunkEnd - 1;
      if (mHalo > 0) {
        while (true) {
          if (data[tmpPtr] != '\n') {
            haloCount++;
            // std::cerr << (char)data[tmpPtr] << std::endl;
            if (haloCount == mHalo) break;
          }
          tmpPtr--;
        }
        // std::cerr << "haloCount: " << haloCount << std::endl;
      } else {
        tmpPtr = chunkEnd;  // nocopy
      }
      std::copy(data + tmpPtr, data + cbufSize, swapBuffer.Pointer());
      bufferSize = cbufSize - tmpPtr;

    } else  // at the end of file
    {
      chunk_->size += r - 1;  // skip the last EOF symbol
      if (usesCrlf) chunk_->size -= 1;

      // only for get seqsinfo
      uint64 chunkEnd = FindCutPos(dataChunk_, data, chunk_->size, mHalo, seqInfos);
      // debug only
      // std::string content((char*)data, chunk_->size);
      // std::cout << "tail content ori: " << data << std::endl;
      // std::cout << "tail content    : " << content << std::endl;
      // end debug

      eof = true;
    }
  } else {
    eof = true;
  }

  return true;
}

uint64 FastaFileReader::FindCutPos(FastaChunk *dataChunk_, uchar *data_, const uint64 size_, const uint64 halo_,
                                   SeqInfos &seqInfos) {
  int count = 0;
  uint64 pos_ = 0;
  uint64 cut_ = 0;       // cut_ point to next '>'
  uint64 lastSeq_ = 0;   //-> the start of last sequences content
  uint64 lastName_ = 0;  //-> the last '>'
  OneSeqInfo seqInfo;
  /*
  [haoz:] dataChunk_ -> start和dataChunk_ -> end 是表示当前这个chunk的开始部分的序列的index和
      最终部分序列的index
  */
  if (data_[0] == '>')  // start with '>'
  {
    dataChunk_->start = this->totalSeqs;
    while (pos_ < size_) {
      if (data_[pos_] == '>') {
        lastName_ = pos_;
        if (FindEol(data_, pos_, size_))  // find name
        {
          ++pos_;
          lastSeq_ = pos_;

          seqInfo.gid = this->totalSeqs;
          seqInfos.push_back(seqInfo);

          this->totalSeqs++;
        } else {
          cut_ = pos_;  // find a cut: incomplete name
          // std::cerr << "cut char: " << (char)data_[cut_] << std::endl;
          break;
        }
      } else {
        ++pos_;
      }
    }

  } else {  // start with {ACGT}
    dataChunk_->start = this->totalSeqs - 1;
    while (pos_ < size_) {
      if (data_[pos_] == '>') {
        lastName_ = pos_;
        if (FindEol(data_, pos_, size_))  // find name
        {
          ++pos_;
          lastSeq_ = pos_;
          seqInfo.gid = this->totalSeqs;
          seqInfos.push_back(seqInfo);
          this->totalSeqs++;
        } else {
          cut_ = pos_;  // find a cut -> '>'
          // std::cout << "cut char: " << (char)data_[cut_] << std::endl;
          break;
        }
      } else {
        ++pos_;
      }
    }
  }

  // no tail cut
  // make sure that cbufSize > name_len + halo
  if (cut_ == 0) {
    uint64 lastSeqLen_ = size_ - lastSeq_;
    if (lastSeqLen_ < halo_) {
      cut_ = lastName_;
      this->totalSeqs--;
    }
  }

  dataChunk_->nseqs = this->totalSeqs - dataChunk_->start;
  dataChunk_->end = this->totalSeqs - 1;

  // if(cut_ != 0) std::cout << "cut: " << cut_ << std::endl;

  return cut_ ? cut_ : size_;
}
bool find_next_seq_start(uchar *data, uint64 size, uint64 &pos_) {
  if (pos_ == size - 1) return false;
  do {
    pos_++;
  } while (pos_ < size && data[pos_] != '>');

  return data[pos_] == '>' ? true : false;
}

// bool FastaFileReader::ReadNextFaChunk(FastaChunk* dataChunk_, SeqInfos& seqInfos, bool &continue_read){
bool FastaFileReader::ReadNextFaChunk(FastaDataChunk *chunk_, SeqInfos &seqInfos, bool &continue_read) {
  if (Eof()) {
    chunk_->size = 0;
    return false;
  }

  // flush the data from previous incomplete chunk
  uchar *data = chunk_->data.Pointer();
  const uint64 cbufSize = chunk_->data.Size();
  chunk_->size = 0;
  int64 toRead = cbufSize - bufferSize;  // buffersize: size left from last chunk

  if (bufferSize > 0) {
    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
    chunk_->size = bufferSize;
    bufferSize = 0;
  }

  // read the next chunk
  int64 r = this->Read(data + chunk_->size, toRead);
  // std::cout << "r is :" << r << std::endl;
  // std::cout << "toRead: " << toRead << std::endl;

  if (r > 0) {
    if (r == toRead)  // somewhere before end
    {
      if (data[0] == '>') {
        uint64 chunkEnd = 0, pos_ = 0;
        pos_ = chunk_->size;
        while (find_next_seq_start(data, cbufSize, pos_)) {
          chunkEnd = pos_;
        }
        if (chunkEnd == 0) {  // means this chunk is a big one, one chunk can not contain all
          continue_read = true;
          chunkEnd = cbufSize;
        } else {
          continue_read = false;
        }
        chunk_->size = chunkEnd - 1;
        if (usesCrlf) chunk_->size -= 1;
        std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
        bufferSize = cbufSize - chunkEnd;
      } else {  // means this is a continue chunk for last sequence
        uint64 chunkEnd = 0, pos_ = 0;
        pos_ = chunk_->size;
        if (find_next_seq_start(data, cbufSize, pos_)) {
          chunkEnd = pos_;
        }
        if (chunkEnd == 0) {  // means this chunk is a big one, one chunk can not contain all
          continue_read = true;
          chunkEnd = cbufSize;
        } else {
          continue_read = false;
        }
        chunk_->size = chunkEnd - 1;
        if (usesCrlf) chunk_->size -= 1;
        std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
        bufferSize = cbufSize - chunkEnd;
      }
    } else {                  // at the end of file
      chunk_->size += r - 1;  // skip the last eof symbol
      if (usesCrlf) chunk_->size -= 1;
      eof = true;
    }
  } else {
    eof = true;
  }

  return true;
}

}  // namespace fa

namespace fq {

int64 count_line(uchar *contenx, int64 read_bytes) {
  int64 count_n = 0;
  for (int i = 0; i < read_bytes; ++i) {
    // printf("%c",contenx[i]);
    if (contenx[i] == '\n') count_n++;
  }
  return count_n;
}

void FastqFileReader::readChunk() {
  FastqDataChunk *part = NULL;

  recordsPool.Acquire(part);
  printf("fastqio: ready to into while\n");

  // while (!errorHandler.IsError() && fileReader.ReadNextChunk(part))
  while (ReadNextChunk(part)) {
    ASSERT(part->size > 0);

    // recordsQueue.Push(numParts, part); //[haoz:] Push:把<numparts,
    // part>组成pair然后放到recordsQueue里面然后notifiy_one
    printf("numParts is %d\n", numParts);
    numParts++;

    recordsPool.Release(part);
    recordsPool.Acquire(part);
  }

  ASSERT(part->size == 0);
  recordsPool.Release(part);  // the last empty part

  // recordsQueue.SetCompleted();
}

FastqDataChunk *FastqFileReader::readNextChunk() {
  FastqDataChunk *part = NULL;
  recordsPool.Acquire(part);
  if (ReadNextChunk(part)) {
    return part;
  } else {
    recordsPool.Release(part);
    return NULL;
  }
}

FastqDataPairChunk *FastqFileReader::readNextPairChunk() {
  FastqDataPairChunk *pair = new FastqDataPairChunk;

  FastqDataChunk *leftPart = NULL;
  recordsPool.Acquire(leftPart);  // song: all in one datapool

  FastqDataChunk *rightPart = NULL;
  recordsPool.Acquire(rightPart);

  int64 left_line_count = 0;
  int64 right_line_count = 0;
  int64 chunkEnd = 0;
  int64 chunkEnd_right = 0;

  //---------read left chunk------------
  if (eof) {
    leftPart->size = 0;
    rightPart->size = 0;
    // return false;
    recordsPool.Release(leftPart);
    recordsPool.Release(rightPart);
    return NULL;
  }

  // flush the data from previous incomplete chunk
  uchar *data = leftPart->data.Pointer();
  const uint64 cbufSize = leftPart->data.Size();
  leftPart->size = 0;
  int64 toRead;
  toRead = cbufSize - bufferSize;
  if (bufferSize > 0) {
    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
    leftPart->size = bufferSize;
    bufferSize = 0;
  }
  int64 r;
  r = Read(data + leftPart->size, toRead);
  if (r > 0) {
    if (r == toRead) {
      chunkEnd =
        cbufSize -  SwapBufferSize; // (1 << 13);  // SwapBufferSize; // SwapBuffersize defined in FastqStream.h as constant value : 1<<13;
      chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
    } else {
      // chunkEnd = r;
      leftPart->size += r - 1;
      if (usesCrlf) leftPart->size -= 1;
      eof = true;
    }
  } else {
    eof = true;
    return NULL;
  }
  //------read left chunk end------//

  //-----------------read right chunk---------------------//
  uchar *data_right = rightPart->data.Pointer();
  const uint64 cbufSize_right = rightPart->data.Size();
  rightPart->size = 0;
  toRead = cbufSize_right - bufferSize2;
  if (bufferSize2 > 0) {
    std::copy(swapBuffer2.Pointer(), swapBuffer2.Pointer() + bufferSize2, data_right);
    rightPart->size = bufferSize2;
    bufferSize2 = 0;
  }
  r = Read2(data_right + rightPart->size, toRead);
  if (r > 0) {
    if (r == toRead) {
      chunkEnd_right =
        cbufSize_right -  SwapBufferSize;
       // (1 << 13);  // SwapBufferSize; // SwapBuffersize defined in FastqStream.h as constant value : 1<<13;
      chunkEnd_right = GetNextRecordPos(data_right, chunkEnd_right, cbufSize_right);
    } else {
      // chunkEnd_right += r;
      rightPart->size += r - 1;
      if (usesCrlf) rightPart->size -= 1;
      eof = true;
    }
  } else {
    eof = true;
    return NULL;
  }
  //--------------read right chunk end---------------------//

  if (!eof) {
    left_line_count = count_line(data, chunkEnd);
    right_line_count = count_line(data_right, chunkEnd_right);
    int64 difference = left_line_count - right_line_count;
    if (difference > 0) {
      // move rightPart difference lines before
      // std::cout << "difference > 0 "<<left_line_count <<" " <<right_line_count << " " << difference <<std::endl;
   //   std::cout << "start: " << chunkEnd_right << std::endl;
      // while(true){
      while (chunkEnd_right < cbufSize) {
        if (data_right[chunkEnd_right] == '\n') {
          difference--;
          if (difference == 0) {
            chunkEnd_right++;
            break;
          }
        }
        chunkEnd_right++;
      }
  //    std::cout << "end: " << chunkEnd_right << std::endl;

    } else if (difference < 0) {
      // move leftPart difference lines before
      // std::cout << "difference < 0 "  <<left_line_count <<" " <<right_line_count << " " << difference << std::endl;
      // while(true){
      while (chunkEnd < cbufSize) {
        if (data[chunkEnd] == '\n') {
          difference++;
          if (difference == 0) {
            chunkEnd++;
            break;
          }
        }
        chunkEnd++;
      }
    }

    leftPart->size = chunkEnd - 1;
    if (usesCrlf) leftPart->size -= 1;

    std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
    bufferSize = cbufSize - chunkEnd;

    rightPart->size = chunkEnd_right - 1;
    if (usesCrlf) rightPart->size -= 1;
    std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer2.Pointer());
    bufferSize2 = cbufSize_right - chunkEnd_right;
  }
  pair->left_part = leftPart;
  pair->right_part = rightPart;
  return pair;
}

FastqDataChunk* FastqFileReader::readNextPairChunkInterleaved(){
  FastqDataChunk *part = NULL;
  recordsPool.Acquire(part);
  uint64 lastChunkPos1;
  uint64 lastChunkPos2;
  //---------read chunk------------
  if (eof) {
    part->size = 0;
    // return false;
    recordsPool.Release(part);
    return NULL;
  }

  // flush the data from previous incomplete chunk
  uchar *data = part->data.Pointer();
  const uint64 cbufSize = part->data.Size();
  part->size = 0;
  int64 toRead;
  toRead = cbufSize - bufferSize;
  if (bufferSize > 0) {
    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
    part->size = bufferSize;
    bufferSize = 0;
  }
  int64 r;
  r = Read(data + part->size, toRead);
  if (r > 0) {
    if (r == toRead) {
      lastChunkPos1 = GetPreviousRecordPos(data, cbufSize - 1, cbufSize);
      lastChunkPos2 = GetPreviousRecordPos(data, lastChunkPos1 - 1, cbufSize);
      
    } else {
      // chunkEnd = r;
      part->size += r - 1;
      if (usesCrlf) part->size -= 1;
      eof = true;
    }
  } else {
    eof = true;
    return NULL;
  }
  //------read chunk end------//
    
  if (!eof) {
    
    part->size = lastChunkPos1 - 1;
    if (usesCrlf) part->size -= 1;

    std::copy(data + lastChunkPos2, data + cbufSize, swapBuffer.Pointer());
    bufferSize = cbufSize - lastChunkPos2;
  }
  return part;

}

// bool IFastqStreamReader::ReadNextChunk(FastqDataChunk* chunk_)
bool FastqFileReader::ReadNextChunk(FastqDataChunk *chunk_) {
  if (Eof()) {
    chunk_->size = 0;
    return false;
  }

  // flush the data from previous incomplete chunk
  uchar *data = chunk_->data.Pointer();
  const uint64 cbufSize = chunk_->data.Size();
  chunk_->size = 0;
  int64 toRead = cbufSize - bufferSize;
  //---------------------------------
  if (bufferSize > 0) {
    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
    chunk_->size = bufferSize;
    bufferSize = 0;
  }

  // read the next chunk
  int64 r = Read(data + chunk_->size, toRead);
  // std::cout << "r is :" << r << std::endl;

  if (r > 0) {
    if (r == toRead)  // somewhere before end
    {
      uint64 chunkEnd = cbufSize - SwapBufferSize;  // Swapbuffersize: 1 << 13
      // std::cout << "chunkend  cbufsize Swapbuffersize: " << chunkEnd <<" "<< cbufSize << " " << SwapBufferSize <<
      // std::endl;
      chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
      chunk_->size = chunkEnd - 1;
      if (usesCrlf) chunk_->size -= 1;

      std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
      bufferSize = cbufSize - chunkEnd;
    } else  // at the end of file
    {
      chunk_->size += r - 1;  // skip the last EOF symbol
      if (usesCrlf) chunk_->size -= 1;

      eof = true;
    }
  } else {
    eof = true;
  }

  return true;
}

uint64 FastqFileReader::GetPreviousRecordPos(uchar *data_, uint64 pos_, const uint64 size_) {
  int offset = 2;

  SkipToSol(data_, pos_, size_);
  if (usesCrlf) {
    offset = 3;
  }
  while (data_[pos_ + offset] != '@') {  //+2
    // std::cout<<"pos_"<<pos_<<std::endl;
    SkipToSol(data_, pos_, size_);
  }
  //标记一下，看是否是质量分
  uint64 pos0 = pos_ + offset;
  SkipToSol(data_, pos_, size_);
  if (data_[pos_ + offset] == '+') {
    //说明上一个@号是质量分
    SkipToSol(data_, pos_, size_);
    SkipToSol(data_, pos_, size_);
    //此时应该是name
    if (data_[pos_ + offset] != '@') {
      std::cout << "core dump is " << data_[pos_ + offset] << std::endl;
      return pos_ + offset;
    } else {
      return pos_ + offset;
    }
  } else {
    return pos0;
  }
}

uint64 FastqFileReader::GetNextRecordPos(uchar *data_, uint64 pos_, const uint64 size_) {
  SkipToEol(data_, pos_, size_);
  ++pos_;

  // find beginning of the next record
  while (data_[pos_] != '@') {
    SkipToEol(data_, pos_, size_);
    ++pos_;
  }
  uint64 pos0 = pos_;

  SkipToEol(data_, pos_, size_);
  ++pos_;

  if (data_[pos_] == '@')  // previous one was a quality field
    return pos_;
  //-----[haoz:] is the following code necessary??-------------//
  SkipToEol(data_, pos_, size_);
  ++pos_;
  if (data_[pos_] != '+') std::cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << std::endl;
  ASSERT(data_[pos_] == '+');  // pos0 was the start of tag
  return pos0;
}

}  // namespace fq

}  // namespace mash
