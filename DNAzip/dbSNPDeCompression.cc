#include "dbSNPDeCompression.h"

vector<bool> _getVecBool(int n, int min_len = 0)
{
	vector<bool> res;
	do
	{
		res.push_back(n%2);
		n/=2;
	} while(n);

	while(res.size() < min_len)
		res.push_back(0);
	reverse(res.begin(), res.end());
	return res;
}

void _prependVec(vector<bool> & dest, vector<bool> src)
{
	dest.insert(dest.begin(), src.begin(), src.end());
}

void _appendVec(vector<bool> & dest, vector<bool> src)
{
	dest.insert(dest.end(), src.begin(), src.end());
}



int _toNum(vector<bool> vec)
{
	int num = 0, f = 1;
	for(int i = vec.size()-1; i >= 0; i--)
	{
		if(vec[i])
			num += f;
		f <<= 1;
	}
	return num;
}

void _printVec(vector<bool> vec, bool new_line = 0)
{
	for(int i = 0; i < vec.size(); i++)
		cout << vec[i];
	if(new_line)
		cout << endl;
}

vector<bool> _fetchVec(vector<bool> & dest, int len)
{
	vector<bool> ans = vector<bool>(dest.begin(), dest.begin() + len);
	dest.erase(dest.begin(), dest.begin() + len);
	return ans;
}

void _decode(vector<bool> & encoded_str, vector<bool> & decoded_str, map< vector<bool>, vector<bool> > encoding)
{
	vector<bool> temp;

	int c = 0;

	for(int i = 0; i < encoded_str.size(); i++)
	{
		temp.push_back(encoded_str[i]);
		if(encoding.find(temp) != encoding.end())
		{
			_appendVec(decoded_str, encoding[temp]);
			temp.clear();
			continue;
		}
	}
	if(temp.size() != 0)
		_appendVec(decoded_str, encoding[temp]);
}

void huffmanDecode(vector<bool> &encoded_str, vector<bool> &decoded_str, int kmer_len = 6)
{
	//cout << "Decompression begins" << endl;
	int dict_size = _toNum(_fetchVec(encoded_str, 10));
	//cout << "Dict Size: " << dict_size << endl;
	vector<bool> new_dict = _fetchVec(encoded_str, dict_size);
	
	map< vector<bool>, vector<bool> > encoding;
	for(int i = 0; i < (1<<kmer_len); i++)
	{
		vector<bool> vec = _getVecBool(i,kmer_len);
		int len = _toNum(_fetchVec(new_dict, 4));
			
		vector<bool> rep = _fetchVec(new_dict, len);
		encoding[rep] = vec;
	}  

	int bits_padding = _toNum(_fetchVec(encoded_str, 3));
	//cout << "Padding: " << bits_padding << endl;

	
	_decode(encoded_str, decoded_str, encoding);
	//printVec(decoded_string,1);
	decoded_str.erase(decoded_str.end() - bits_padding, decoded_str.end());
	
	//cout << decoded_str.size() << endl;
}



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
		
		/* Decoding the Huffman Encoding of the Bitvector */
		vector<bool> encoded_str;
		for(int i = 0; i < bitMapSz; i++)
			encoded_str.push_back(bitMap[i]);
		
		vector<bool> decoded_str;
		huffmanDecode(encoded_str, decoded_str);
		cout << "Decoded BitVector From " << bitMapSz << " to " << decoded_str.size() << endl;
		
		delete bitMap;
		bitMapSz = decoded_str.size();
		bitMap = new bool[bitMapSz];
		for(int i = 0; i < bitMapSz; i++)
			bitMap[i] = decoded_str[i];
		
		/* Decoding ends */	
				
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
		std::map<int, string>* posGen = new map<int, string>();	 //pos & genLetter	
		vector<int>* positions = new vector<int>(); //vector of positions								
			
		//a bitMap of common DELs				
		unsigned bitExists = readBitVINT(srcBf);
		bool *bitMap;
		unsigned bitMapSz;

		if(bitExists) {
			bitMapSz = readBitVINT( srcBf );
			bitMap = new bool[bitMapSz];
			readBitArrays( srcBf, bitMap, bitMapSz );
		
		
			/* Decoding the Huffman Encoding of the Bitvector */
			vector<bool> encoded_str;
			for(int i = 0; i < bitMapSz; i++)
				encoded_str.push_back(bitMap[i]);
			
			vector<bool> decoded_str;
			huffmanDecode(encoded_str, decoded_str, 7);
			cout << "Decoded BitVector From " << bitMapSz << " to " << decoded_str.size() << endl;
			
			delete bitMap;
			bitMapSz = decoded_str.size();
			bitMap = new bool[bitMapSz];
			for(int i = 0; i < bitMapSz; i++)
				bitMap[i] = decoded_str[i];
		
			/* Decoding ends */	
		}		
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

		if(bitExists)
		{
			int bitCnt = 0;
			while( getline( dbSNPStream, dbSNPLine ) )
			{
				if( dbSNPLine.find("deletion") == string::npos )
				continue;
	
				string observed;
				
				strtok( (char *)dbSNPLine.c_str(), "\t" ); //ID					
				strtok( NULL, "\t" ); //chrNum					
				strtok( NULL, "\t" ); //startPos
				int dbSNPPos = atoi( strtok(NULL, "\t") ); //endPos
				strtok( NULL, "\t"); //name
				strtok( NULL, "\t"); //score
				strtok( NULL, "\t"); //strand
				observed = strtok( NULL, "\t"); //refNCBI
				
				if( bitMap[bitCnt++] )
				{	
					(*posGen)[dbSNPPos] = observed;			    
				}
		}
		dbSNPStream.close();		
		delete bitMap;			
		}
		//
		//rare SNPs
		//		
		
		unsigned rareSNPCnt = readBitVINT( srcBf );			
		cout << "Received Count: "<<rareSNPCnt << endl;
		unsigned prevPos = 0;
		for( unsigned i = 0; i < rareSNPCnt; i++ )
		{
			unsigned diff = readBitVINT( srcBf );
			positions->push_back( prevPos + diff );
			prevPos = prevPos + diff;				
		}
		srcBf.ByteAlign();
						
		vector<int>::iterator posIter = positions->begin();
		while( posIter != positions->end() )
		{
				unsigned length = readBitVINT (srcBf);
				//cout << refGen[0] << " " << *posIter << " "<< length << "\n";
			    (*posGen)[*posIter] = string( refGen + *posIter - length, length);
				posIter++;
		}
			
		map<int, string>::iterator it;			
		for ( it=posGen->begin(); it != posGen->end(); it++ )
		{
			destStream <<DELETION<<","<<chromosomeID <<","<< it->first <<","<< it->second <<"/"<< string(it->second.length(), '-') << endl;				
		}	
						
		delete posGen;
		delete positions;			
		
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

