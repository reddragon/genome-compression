#ifndef WATSONFILEGENERATOR_H_
#define WATSONFILEGENERATOR_H_

#include <fstream>
#include <iostream>
#include "fileGenerator.h"
#include <cstring>

class watsonFileGenerator: public fileGenerator
{
 protected:
	 string SNPFile;
	 string indelFile;
	 string destFile;
 public:
	 watsonFileGenerator( const string& SNPFile,
			              const string& indelFile, 
			              const string& destFile );
     //This forces every derived class to have 
     //its own generate function
     void generate();      
};
#endif /*WATSONFILEGENERATOR_H_*/
