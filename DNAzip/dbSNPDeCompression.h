#ifndef DBSNPNEWDECOMPRESSION_H_
#define DBSNPNEWDECOMPRESSION_H_

#define GENLENGTH 300000000

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include "stdlib.h"
#include "output.h"
#include "input.h"
#include "huffman.h"
#include "fileGenerator.h"
#include <cstring>
//#include "bithuffman.h"

using namespace std;



class dbSNPDeCompression{
private:
	string dbSNPDir;
	string dataDir;
	const signed kMer;
	Hufftree<string, int>* hufftree;
	char* refGen; 
	
	void deCompressSNPs( bit_file_c& srcBf, 
		ofstream& destStream, 
		const string& chromosomeID  );
	void deCompressINs( bit_file_c& srcBf, 
		ofstream& destStream, 
		const string& chromosomeID );
	void deCompressDELs( bit_file_c& srcBf, 
		ofstream& destStream, 
		const string& chromosomeID );
public:
	dbSNPDeCompression( const string& inFreqFile = "files//JWB-in-posFreq4.txt", 
			const string& dbSNPDir = "dbSNP//",
			const string& dataDir = "chr//", 
			unsigned kMer = 4 );
	void deCompress( const string& srcFile, const string& destFile );
	~dbSNPDeCompression();
	
};

#endif /*DBSNPNEWDECOMPRESSION_H_*/
