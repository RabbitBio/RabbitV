#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace std;

static constexpr char reverseArray[8] = {'T', 'G', 'A', 'C', 'N', 'N', 'N', 'N'};
static constexpr uint8_t reverseMask = 0x0e;

class Sequence{
public:
    Sequence();
    Sequence(string seq);
    void print();
    int length();
    Sequence reverseComplement();

    Sequence operator~();

    static bool test();

public:
    string mStr;
};

#endif
