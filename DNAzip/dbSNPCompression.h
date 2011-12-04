#ifndef DBSNPNEWCOMPRESSION_H_
#define DBSNPNEWCOMPRESSION_H_

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <cassert>
#include "stdlib.h"
#include "bitfile.h"
#include "output.h"
#include "input.h"
#include "huffman.h"
#include "fileGenerator.h"
#include <cstring>
#include "bithuffman.h"
#include "rle.h"

using namespace std;

class dbSNPCompression{
private:		
	string dbSNPDir;
	const unsigned kMer;
	Hufftree<string, int>* hufftree;
	
	void compressSNPs( const string& srcFile, bit_file_c& destBf );
	void compressINs( const string& srcFile,  bit_file_c& destBf );
	void compressDELs( const string& srcFile, bit_file_c& destBf );
	
public:
	dbSNPCompression( const string& inFreqFile = "files//WB-in-posFreq4.txt", 
			const string& dbSNPDir = "dbSNP//", 
			unsigned kMer = 4 );
	void compress( const string& srcFile, const string& destFile );
	~dbSNPCompression();	
		
};


#endif /*DBSNPNEWCOMPRESSION_H_*/
