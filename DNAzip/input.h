#ifndef INPUT_H_
#define INPUT_H_

#include <vector>
#include "bitfile.h"

using namespace std;

//read a vector of bits
void readBits( bit_file_c& bf, unsigned count, vector<bool>* bits );
//read an array of bits
void readBitArrays( bit_file_c& bf, bool* bits, unsigned cnt );
//read a vector of gene letters
void readBitGens( bit_file_c& bf, unsigned count, vector<char>* gens );
//read a VINT 
unsigned readBitVINT( bit_file_c& bf );
//read a string
string readString( bit_file_c& bf );


#endif /*INPUT_H_*/
