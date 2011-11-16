#ifndef KMERFREQGENERATOR_H_
#define KMERFREQGENERATOR_H_

#include <fstream>
#include <iostream>
#include <map>
#include "fileGenerator.h"
#include "output.h"
#include <cstring>

class kMerFreqGenerator: public fileGenerator
{
protected:
		int kMer;
		string srcFile;
		string destFile;
 public:
	 kMerFreqGenerator( int kMer, const string& srcFile, const string& destFile );
     //This forces every derived class to have 
     //its own generate function
     void generate();      
};

#endif /*KMERFREQGENERATOR_H_*/
