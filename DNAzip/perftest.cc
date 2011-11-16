#include <map>
#include "input.h"
#include "output.h"
//#include "dbSNPCompression.h"
//#include "dbSNPDeCompression.h"
#include "dbSNPCompression.h"
#include "dbSNPDeCompression.h"
#include "watsonFileGenerator.h"
#include "kMerFreqGenerator.h"


void testKMerFreq( int freq, const string& originalFile, const string& destionationFile )
{
	kMerFreqGenerator gen( freq, originalFile, destionationFile );
	gen.generate();	
	cout << "testKMerFreq succeed" <<endl;	
}

void testdbCompression( const string& kMerFreqFile, const string& originalFile, const string& destionationFile )
{
	dbSNPCompression dbCP(kMerFreqFile);	
	dbCP.compress( originalFile, destionationFile );
	cout << "testdbCompression succeed" <<endl;
}

void testdbDeCompression( const string& kMerFreqFile, const string& originalFile, const string& destionationFile )
{
	dbSNPDeCompression dbDeCP(kMerFreqFile);	
	dbDeCP.deCompress( originalFile, destionationFile );
	cout << "testdbDeCompression succeed" <<endl;
}

void testWatsonFileGenerator( const string& SNPFile,
        const string& indelFile, 
        const string& destFile )
{
	watsonFileGenerator gen( SNPFile, indelFile, destFile );
	gen.generate();	
	cout << "testWatsonFileGenerator succeed" <<endl;
}

int main()
{
	//generate a data file for James Watson Genome Sequence
	testWatsonFileGenerator(  "files//JWB-snps-submission.txt", 
			"files//JWB-indels-submission", "files//JWB-unified-file.txt" );
	//generate the frequencies of all the KMers
	testKMerFreq(4, "files//JWB-unified-file.txt", "files//JWB-in-posFreq4.txt");	
	//do the compression	
	testdbCompression( "files//JWB-in-posFreq4.txt", "files//JWB-unified-file.txt", "files//JWB-unified-compression.txt" );
	//do the deCompression
	testdbDeCompression( "files//JWB-in-posFreq4.txt", "files//JWB-unified-compression.txt", "files//JWB-unified-DeCompression.txt" );	
} 
  
		

