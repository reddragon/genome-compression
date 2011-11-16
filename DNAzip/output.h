#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <vector>
#include "bitfile.h"

using namespace std;

//write a vector of bits
void writeBits( bit_file_c& bf, vector<bool>* bits );
//write an array of bits
void writeBitArrays( bit_file_c& bf, bool* bits, unsigned cnt );
//write a vector of gene letters
void writeBitGens( bit_file_c& bf, vector<char> * gens );
//write a VINT 
void writeBitVINT( bit_file_c& bf, unsigned num );
//write a string
void writeString( bit_file_c& bf, string& str );


#endif /*OUTPUT_H_*/
