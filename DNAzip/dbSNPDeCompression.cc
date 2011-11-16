#include "dbSNPDeCompression.h"

dbSNPDeCompression::dbSNPDeCompression( const string& inFreqFile, 
		const string& dbSNPDir, 
		const string& dataDir, 
		unsigned kMer):dbSNPDir(dbSNPDir), dataDir(dataDir), kMer(kMer)
{
	std::map<string, int>* frequencies = new map<string, int>();
		
	bit_file_c bf; 
	bf.Open( inFreqFile.c_str(), BF_READ );
		
	unsigned num = readBitVINT( bf );  	
	for( unsigned i = 0; i < num; i++ )
	{
		string gens = readString( bf );
		int freq = readBitVINT( bf );
		(*frequencies)[ gens ] = freq;
	}
	bf.Close();
				
	hufftree = new Hufftree<string, int>(frequencies->begin(), frequencies->end());
	delete frequencies;
	
	refGen = new char[GENLENGTH];
}

dbSNPDeCompression::~dbSNPDeCompression()
{
	delete hufftree;
	delete refGen;
}

void dbSNPDeCompression::deCompressSNPs( bit_file_c& srcBf, 
		ofstream& destStream, 
		const string& chromosomeID )
{				
		
		std::map<int, char>* posGen = new map<int, char>();	 //pos & genLetter	
		vector<char>* gens = new vector<char>();  //vector of genLetters
		vector<int>* positions = new vector<int>(); //vector of positions								
			
		//a bitMap of common SNPs				
		
		unsigned bitMapSz = readBitVINT( srcBf );
		bool* bitMap = new bool[bitMapSz];
		readBitArrays( srcBf, bitMap, bitMapSz );								
					
		string dbSNPLine;
		ifstream dbSNPStream( (dbSNPDir+chromosomeID+".txt").c_str() );
		getline( dbSNPStream, dbSNPLine ); //schema of dbSNP
		
		string dataLine;
		ifstream dataStream ( (dataDir+chromosomeID+".fa").c_str() );		
		getline( dataStream, dataLine );   //original data file				
		int genLetterCnt = 0;				
		while( getline(dataStream, dataLine) )
		{			
			for( unsigned j = 0; j < dataLine.size(); j++ )
			{
				refGen[genLetterCnt++] = dataLine[j];
			}			
		}
		dataStream.close();
		
		int bitCnt = 0;
		while( getline( dbSNPStream, dbSNPLine ) )
		{
			if( dbSNPLine.find("single") == string::npos )
			continue;
			
			char* subToken;
			strtok( (char *)dbSNPLine.c_str(), "\t" ); //ID					
			strtok( NULL, "\t" ); //chrNum					
			strtok( NULL, "\t" ); //startPos
			int dbSNPPos = atoi( strtok(NULL, "\t") ); //endPos
			strtok( NULL, "\t"); //name
			strtok( NULL, "\t"); //score
			strtok( NULL, "\t"); //strand
			strtok( NULL, "\t"); //refNCBI
			strtok( NULL, "\t"); //refUCSC
			subToken = strtok( NULL, "\t"); //observed
			
			if( bitMap[bitCnt++] )
			{	
				if( dbSNPPos <= genLetterCnt )
				{
					if(  refGen[dbSNPPos-1] == subToken[0]  )
					{
						(*posGen)[dbSNPPos] = subToken[2];			    
					}
					else 
					{
						(*posGen)[dbSNPPos] = subToken[0];			    
					}
				}
			}
		}
		dbSNPStream.close();		
		delete bitMap;			
		//
		//rare SNPs
		//		
		
		unsigned rareSNPCnt = readBitVINT( srcBf );			
		unsigned prevPos = 0;
		for( unsigned i = 0; i < rareSNPCnt; i++ )
		{
			unsigned diff = readBitVINT( srcBf );
			positions->push_back( prevPos + diff );
			prevPos = prevPos + diff;				
		}
		readBitGens( srcBf, rareSNPCnt, gens );
		srcBf.ByteAlign();
						
		vector<int>::iterator posIter = positions->begin();
		vector<char>::iterator genIter = gens->begin();
		while( posIter != positions->end() )
		{
			    (*posGen)[*posIter] = *genIter;					
				posIter++;
				genIter++;
		}
			
		map<int, char>::iterator it;			
		for ( it=posGen->begin(); it != posGen->end(); it++ )
		{
			destStream <<SNP<<","<<chromosomeID <<","<< it->first <<","<< refGen[it->first-1] <<"/"<< it->second<<endl;				
		}	
						
		delete posGen;
		delete gens;
		delete positions;			
		
		cout << "deCompressSNPs" <<chromosomeID <<endl;
}

void dbSNPDeCompression::deCompressINs( bit_file_c& srcBf, 
		ofstream& destStream, 
		const string& chromosomeID )
{
	vector<unsigned>* positions = new vector<unsigned>(); // a vector of positions
	vector<unsigned>* lengths = new vector<unsigned>();   // a vector of lengths
	vector<string>* longGens = new vector<string>();      // a vector of longGens	
	
	//number of insertions
	unsigned inCnt = readBitVINT( srcBf );
	//a vector of positions
	unsigned prevPos = 0;
	for( unsigned i = 0; i < inCnt; i++ )
	{
		prevPos += readBitVINT( srcBf );			
		positions->push_back( prevPos );			
	}
	//a vector of lengths	
	for( unsigned i = 0; i < inCnt; i++ )
	{
		unsigned length = readBitVINT( srcBf );
		lengths->push_back( length );		
	}
	//number of bits for long gens
	vector<bool> genBits;
	unsigned bitCnt = readBitVINT( srcBf );		
	readBits( srcBf, bitCnt, &genBits );		
	hufftree->decode( genBits, std::back_inserter(*longGens) );	
	
	vector<unsigned>::iterator posIter = positions->begin();
	vector<unsigned>::iterator lenIter = lengths->begin();
	vector<string>::iterator longGenIter = longGens->begin();
	
	while( posIter != positions->end() )
	{
		unsigned pos = *posIter;
		unsigned len = *lenIter;
		string gens = "";
		for( unsigned i = 0; i < (len/kMer); i++ )
		{			
			gens += (*longGenIter);
			longGenIter++;
		}
		if( len%kMer > 0 )
		{
			gens += (*longGenIter);
			longGenIter++;
		}
		
		char* refStr = new char[len+1];
		for(unsigned j = 0; j < len; j++)
		{
			refStr[j] = '-';
		}
		refStr[len] = '\0';
		destStream<< INSERTION <<","<<chromosomeID << "," << pos << ","<< refStr<<"/"<<gens <<endl;
			
		posIter++;
		lenIter++;
	}	
		
	srcBf.ByteAlign();
		
	delete positions;
	delete lengths;
	delete longGens;	
		
}

void dbSNPDeCompression::deCompressDELs( bit_file_c& srcBf, 
		ofstream& destStream, 
		const string& chromosomeID )
{
	string dataLine;
	ifstream dataStream ( (dataDir+chromosomeID+".fa").c_str() );		
	getline( dataStream, dataLine );   //original data file				
	int genLetterCnt = 0;	
	while( getline(dataStream, dataLine) )
	{		
		for( unsigned j = 0; j < dataLine.size(); j++ )
	    {
			refGen[genLetterCnt++] = dataLine[j];
		}			
	}
	dataStream.close();
	
	vector<unsigned>* positions = new vector<unsigned>(); // positions
	vector<unsigned>* lengths = new vector<unsigned>();   // lengths				
	
	//number of deletions
	unsigned delCnt = readBitVINT( srcBf );
	unsigned prevPos = 0;
	for( unsigned i = 0; i < delCnt; i++ )
	{
		prevPos += readBitVINT( srcBf );			
		positions->push_back( prevPos );			
	}			
	for( unsigned i = 0; i < delCnt; i++ )
	{
		unsigned length = readBitVINT( srcBf );
		lengths->push_back( length );				
	}							
		
	vector<unsigned>::iterator posIter = positions->begin();
	vector<unsigned>::iterator lenIter = lengths->begin();
		
	while( posIter != positions->end() )
	{
		char* refStr = new char[(*lenIter)+1];
		char* altStr = new char[(*lenIter)+1];
		for( unsigned j = 0; j < (*lenIter); j++ )
		{
			refStr[j] = refGen[(*posIter)-1+j];
			altStr[j] = '-';
		}
		refStr[(*lenIter)] = '\0';
		altStr[(*lenIter)] = '\0';
			
		destStream << DELETION << "," << chromosomeID << "," << (*posIter) << ","<< refStr <<"/"<<altStr<<endl;
		
		delete refStr;
		delete altStr;
		
		posIter++;
		lenIter++;
	}	
		
	srcBf.ByteAlign();
			
	delete positions;
	delete lengths;	  	
	
	cout << "deCompressDELs" <<chromosomeID <<endl;
}

void dbSNPDeCompression::deCompress( const string& srcFile, 
		const string& destFile )
{
	string chromosomeID;
	bit_file_c srcBf;
	srcBf.Open( srcFile.c_str(), BF_READ );
	
	ofstream destStream( destFile.c_str() );
	
	while( (chromosomeID = readString ( srcBf )) != "" )
	{
		unsigned opType = readBitVINT ( srcBf );
		if( opType == SNP )
			deCompressSNPs( srcBf, destStream, chromosomeID );
		else if( opType == INSERTION )
			deCompressINs( srcBf, destStream, chromosomeID );
		else
			deCompressDELs( srcBf, destStream, chromosomeID );			
	}
	
	srcBf.Close();
	destStream.close();
}

