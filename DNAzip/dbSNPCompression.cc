#include "dbSNPCompression.h"

dbSNPCompression::dbSNPCompression( const string& inFreqFile, 
		const string& dbSNPDir, 
		unsigned kMer ): dbSNPDir(dbSNPDir), kMer(kMer)
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
	   
}

dbSNPCompression::~dbSNPCompression()
{	
	delete hufftree;
}

void dbSNPCompression::compressSNPs( 
		const string &srcFile, 
		bit_file_c& destBf )
{
	    ifstream srcStream( srcFile.c_str() );
		string line;			
																
		std::map<int, char>* refPosGen = new map<int, char>();  //pos & genLetter in ref
		std::map<int, char>* altPosGen = new map<int, char>();  //pos & genLetter in alt
		
		vector<char>* gens = new vector<char>();   //vector of genLetters
		vector<int>* positions = new vector<int>();  //vector of positions		
		
		string prevChrID = "";							
		unsigned prevPos = 0;

		while( getline( srcStream, line ) )
		{				 			
			//Tokenize the line						
			char* token;
			token = strtok ( (char *)line.c_str(), "," ); //type
			int type = atoi( token );
			if( type != SNP )
				continue;
			
			token = strtok ( NULL, "," ); //chrNum
			string newChrID = token;			
			token = strtok ( NULL, "," ); //position
			int newPos = atoi(token);			
			token = strtok ( NULL, "," ); //genLetter		
			char newGen = token[2];
			char oldGen = token[0];
								
			if( newChrID.compare(prevChrID) != 0 )
			{
				if(  refPosGen->size() > 0 )
				{								
					//create a bitVector of dbSNP
					string subLine;
					ifstream dbSNPCntStream( (dbSNPDir+prevChrID+".txt").c_str() );
					getline( dbSNPCntStream, subLine ); //schema of dbSNP
					unsigned bitCnt = 0;
					while( getline( dbSNPCntStream, subLine ) )
					{
						if( subLine.find("single") == string::npos )
							continue;
						bitCnt++;
					}									
					bool* bitMap = new bool[bitCnt];								
					for( unsigned i = 0; i < bitCnt; i++ )					
						bitMap[i] = false;					
					dbSNPCntStream.close();
				
					ifstream dbSNPStream( (dbSNPDir+prevChrID+".txt").c_str() );
					getline( dbSNPStream, subLine ); //schema of dbSNP				
					bitCnt = 0;
					while( getline( dbSNPStream, subLine ) )
					{
						if( subLine.find("single") == string::npos )
							continue;
						char* subToken;
						strtok( (char *)subLine.c_str(), "\t" ); 	//ID					
						strtok( NULL, "\t" ); 	//chrNum					
						strtok( NULL, "\t" ); 	//startPos
						int dbSNPPos = atoi( strtok(NULL, "\t") );  //endPos					
						map<int, char>::iterator refIt;
						map<int, char>::iterator altIt;
						if( (refIt = refPosGen->find(dbSNPPos)) !=  refPosGen->end() )
						{
							altIt = altPosGen->find( dbSNPPos );
						
							strtok( NULL, "\t"); //name
							strtok( NULL, "\t"); //score
							strtok( NULL, "\t"); //strand
							strtok( NULL, "\t"); //refNCBI
							strtok( NULL, "\t"); //refUCSC
							subToken = strtok( NULL, "\t"); //observed
							
							if( ( subToken[0] == refIt->second && subToken[2] == altIt->second )
							|| ( subToken[2] == refIt->second && subToken[0] == altIt->second ) )
							{
								bitMap[bitCnt] = true;							
								refPosGen->erase( refIt );
								altPosGen->erase( altIt );								
							}												
							else
							{
								bitMap[bitCnt] = false;							
							}
						}
						bitCnt++;						 					 
					}
					dbSNPStream.close();
					//chrID									
					writeString( destBf, prevChrID);
					//operation type
					writeBitVINT( destBf, SNP );
					//size of bitmap
					writeBitVINT( destBf, bitCnt );
					//bitMap
					writeBitArrays( destBf, bitMap, bitCnt );											
				
					delete bitMap;				
				
					//create a vector for other SNPs
					map<int, char>::iterator it;
					int pos = 0; 
					for ( it=altPosGen->begin() ; it != altPosGen->end(); it++ )
					{
						positions->push_back( (it->first - pos) );
						gens->push_back( it->second );
						pos = it->first;
					}												
					//size of other positions				
					writeBitVINT( destBf, positions->size() );
					//positions for other SNPs
					vector<int>::iterator posIter = positions->begin();				
					while( posIter != positions->end() )
					{					
						writeBitVINT( destBf, (unsigned)(*posIter) );
						posIter++;
					}						
					//gens
					writeBitGens(destBf, gens);
					destBf.ByteAlign();															
												
					positions->clear();
					gens->clear();
					refPosGen->clear();
					altPosGen->clear();
				}
				
				prevChrID = newChrID;
				prevPos = 0;							
			}												
			
			// store the positions & gens into refPosGen & altPosGen
			if( newPos > (int)prevPos )
			{
			   (*refPosGen)[ newPos ] = oldGen;
			   (*altPosGen)[ newPos ] = newGen;
				prevPos = newPos;																						
			}						
		}		
		
		//compress the last chromosome
		//create a bitVector of dbSNP
		string subLine;
		ifstream dbSNPCntStream( (dbSNPDir+prevChrID+".txt").c_str() );
		getline( dbSNPCntStream, subLine ); //schema of dbSNP
		unsigned bitCnt = 0;
		while( getline( dbSNPCntStream, subLine ) )
		{
			if( subLine.find("single") == string::npos )
				continue;
			bitCnt++;
		}
		bool* bitMap = new bool[bitCnt];						
		for( unsigned i = 0; i < bitCnt; i++ )
		{
			bitMap[i] = false;
		}			
		
		ifstream dbSNPStream( (dbSNPDir+prevChrID+".txt").c_str() );
		getline( dbSNPStream, subLine ); //schema of dbSNP								
		bitCnt = 0;
		while( getline( dbSNPStream, subLine ) )
		{
			if( subLine.find("single") == string::npos )
				continue;
			char* subToken;
			strtok( (char *)subLine.c_str(), "\t" ); //ID
			strtok( NULL, "\t" ); //chrNum
			strtok( NULL, "\t" ); //startPos
			int dbSNPPos = atoi( strtok(NULL, "\t") ); //endPos					
			map<int, char>::iterator refIt;
			map<int, char>::iterator altIt;					
			if( (refIt = refPosGen->find(dbSNPPos)) != refPosGen->end() )
			{
				altIt = altPosGen->find( dbSNPPos );
				
				strtok( NULL, "\t"); //name
		    	strtok( NULL, "\t"); //score
				strtok( NULL, "\t"); //strand
				strtok( NULL, "\t"); //refNCBI
				strtok( NULL, "\t"); //refUCSC
				subToken = strtok( NULL, "\t"); //observed
				if( ( subToken[0] == refIt->second && subToken[2] == altIt->second )
				|| ( subToken[2] == refIt->second && subToken[0] == altIt->second ) )
				{
					bitMap[bitCnt] = true;							
					refPosGen->erase( refIt );
					altPosGen->erase( altIt );					
				}		
				else
				{
					bitMap[bitCnt] = false;							
				}
			}
			bitCnt++;						 					 
		}
		dbSNPStream.close();
		
		//chrID						
		writeString( destBf, prevChrID);
		//operation type
		writeBitVINT( destBf, SNP );
		//size of bitmap
		writeBitVINT( destBf, bitCnt );
		//bitMap
		writeBitArrays( destBf, bitMap, bitCnt );		
		delete bitMap;		
		
		//create a vector for other SNPs				
		map<int, char>::iterator it;
		int pos = 0; 
		for ( it = altPosGen->begin() ; it != altPosGen->end(); it++ )
		{
		    positions->push_back( (it->first - pos) );
			gens->push_back( it->second );
			pos = it->first;
		}												
		//size of other SNPs
		writeBitVINT( destBf, positions->size() );
		//positions for other SNPs
		vector<int>::iterator posIter = positions->begin();				
		while( posIter != positions->end() )
		{					
			writeBitVINT( destBf, (unsigned)(*posIter) );
			posIter++;
		}						
		//gens
		writeBitGens(destBf, gens);	
		destBf.ByteAlign();							

		delete refPosGen;
		delete altPosGen;
		delete gens;
		delete positions;
		
		srcStream.close();				
}

void dbSNPCompression::compressINs( 
		const string &srcFile, 
		bit_file_c& destBf )
{
	ifstream srcStream( srcFile.c_str() );
	string line;							
									
	vector<unsigned>* positions = new vector<unsigned>(); // vector for positions
	vector<unsigned>* lengths = new vector<unsigned>(); // vector for lengths
	vector<string>* longGens = new vector<string>(); // vector for longGens whose length is larger than kMer			
			
	string prevChrID = "";								
	signed oldPos = 0;

	while( getline( srcStream, line ) )
	{				 			
		//Tokenize the line						
		char* token;
		token = strtok ( (char *)line.c_str(), "," ); //type
		int type = atoi( token );
		if( type != INSERTION )
			continue;
		token = strtok ( NULL, "," ); //chromsome ID
		string newChrID = token; 
				
		token = strtok ( NULL, "," ); //pos		
		signed newPos = atoi(token);  		
						
		token = strtok ( NULL, "," ); //genLetters
		string genToken = token;
													
		if( newChrID.compare(prevChrID) != 0 )
		{
			if(  positions->size() > 0 )
			{
				//chromsome ID
				writeString( destBf, prevChrID );
				//operation type
				writeBitVINT( destBf, INSERTION );
				//total number of insertions
				writeBitVINT( destBf, positions->size() );				
				//positions
				vector<unsigned>::iterator posIter = positions->begin();
				unsigned prevPos = 0;
				while( posIter != positions->end() )
				{
					unsigned pos = *posIter;
					writeBitVINT( destBf, (unsigned)(pos-prevPos) );
					prevPos = pos;
					posIter++;
				}				
				//lengths			
				vector<unsigned>::iterator lenIter = lengths->begin();			 
				while( lenIter != lengths->end() )
				{
					unsigned len = *lenIter;
					writeBitVINT( destBf, len );			 		
					lenIter++;
				}					
				//encoded long gens
				vector<bool> encodedLongGens = hufftree->encode(
				longGens->begin(), longGens->end() ); 				
				writeBitVINT( destBf, encodedLongGens.size() );
				writeBits( destBf, &encodedLongGens );											
				destBf.ByteAlign();
				//clear
				positions->clear();
				lengths->clear();
				longGens->clear();											
			 }
			 prevChrID = newChrID;
			 oldPos = 0;			 
		}			
		if( newPos > oldPos )
		{
			positions->push_back( newPos );			
			signed strLen = genToken.length()/2;
			lengths->push_back( strLen );
			string strGens = genToken.substr( strLen+1, strLen );
			signed numOfLongs = strGens.length()/kMer;
			for( signed i = 0; i < numOfLongs; i++ )
			{				
				longGens->push_back( strGens.substr( i*kMer, kMer ) );
			}
			if( strGens.length()%kMer > 0 )					
			{				
				longGens->push_back( strGens.substr( numOfLongs*kMer, strGens.length()%kMer ) );
			}
			oldPos = newPos;	
		}
	}
		
	if(  positions->size() > 0 )
	{
		//chromsome ID
		writeString( destBf, prevChrID );
		//operation type
		writeBitVINT( destBf, INSERTION );
		//total number of insertions
		writeBitVINT( destBf, positions->size() );
		//positions
		vector<unsigned>::iterator posIter = positions->begin();
		unsigned prevPos = 0;
		while( posIter != positions->end() )
		{
			unsigned pos = *posIter;
			writeBitVINT( destBf, (unsigned)(pos-prevPos) );
			prevPos = pos;
			posIter++;
		}
		//lengths			
		vector<unsigned>::iterator lenIter = lengths->begin();			 
		while( lenIter != lengths->end() )
		{
			unsigned len = *lenIter;
			writeBitVINT( destBf, len );			 		
			lenIter++;
		}						
		//encoded long gens
		vector<bool> encodedLongGens = hufftree->encode(
		longGens->begin(), longGens->end()); 
		writeBitVINT( destBf, encodedLongGens.size() );
		writeBits( destBf, &encodedLongGens);		
		destBf.ByteAlign();
				
	}
	//clear
	delete positions;
	delete lengths;
	delete longGens;	
					
	srcStream.close();			
}

void dbSNPCompression::compressDELs( 
		const string &srcFile, 
		bit_file_c& destBf )
{
	ifstream srcStream(  srcFile.c_str() );
	string line;								
		
	vector<unsigned>* positions = new vector<unsigned>(); //positions
	vector<unsigned>* lengths = new vector<unsigned>(); //lengths	
		
	string prevChrID = "";	
	signed oldPos = 0;
	
	while( getline( srcStream, line ) )
	{						
		//Tokenize the line						
		char* token;
		token = strtok ( (char *)line.c_str(), "," ); //type
		int type = atoi( token );
		if( type != DELETION )
			continue;
		token = strtok ( NULL, "," ); //chromsome ID
		string newChrID = token; 
						
		token = strtok ( NULL, "," ); //pos		
		signed newPos = atoi(token);  		
								
		token = strtok ( NULL, "," ); //genLetters
		string genToken = token;
					
		if( newChrID.compare(prevChrID) != 0 )
		{
			if( positions->size() > 0 )
			{				
				//chromsome ID
				writeString( destBf, prevChrID );
				//operation type
				writeBitVINT( destBf, DELETION );
				//size of positions
				writeBitVINT( destBf, positions->size() );
				//positions
				vector<unsigned>::iterator posIter = positions->begin();
				unsigned prevPos = 0;
				while( posIter != positions->end() )
				{
					unsigned newPos = *posIter;
					writeBitVINT( destBf, (unsigned)(newPos-prevPos) );
					prevPos = newPos;
					posIter++;
				}
				//lengths
				vector<unsigned>::iterator lenIter = lengths->begin();			 
				while( lenIter != lengths->end() )
				{
				 	unsigned len = *lenIter;
				 	writeBitVINT( destBf, len );			 		
					lenIter++;
				}
				destBf.ByteAlign();
				
				positions->clear();
				lengths->clear();
			}	
			
			prevChrID = newChrID;
			oldPos = 0;						
		}
		
		if( newPos > oldPos )
		{
			unsigned genLen = genToken.length()/2; //genLength					
			positions->push_back( newPos );
			lengths->push_back( genLen );					
			oldPos = newPos;
		 }						
	}	 
	
	if( positions->size() > 0 )
	{		
		//chromsome ID
		writeString( destBf, prevChrID );
		//operation type
		writeBitVINT( destBf, DELETION );
		//size of positions
		writeBitVINT( destBf, positions->size() );
		//positions
		vector<unsigned>::iterator posIter = positions->begin();
		unsigned prevPos = 0;
		while( posIter != positions->end() )
		{
			unsigned newPos = *posIter;
			writeBitVINT( destBf, (unsigned)(newPos-prevPos) );
			prevPos = newPos;
			posIter++;
		}
		//lengths
		vector<unsigned>::iterator lenIter = lengths->begin();			 
		while( lenIter != lengths->end() )
		{
			 unsigned len = *lenIter;
			 writeBitVINT( destBf, len );			 		
		     lenIter++;
		}
		destBf.ByteAlign();									
	}	
						
	delete positions;
	delete lengths;	 
			
	srcStream.close();		
}

void dbSNPCompression::compress( const string& srcFile, 
		const string& destFile )
{
	bit_file_c destBf;
	destBf.Open( destFile.c_str(), BF_WRITE );
	//SNPs
	compressSNPs( srcFile, destBf );
	cout << "compress SNPs" <<endl;	
	//deletions
	compressDELs( srcFile, destBf );
	cout << "compress DELs" <<endl;
	//insertions
	compressINs( srcFile, destBf );
	cout << "compress INs" <<endl;
	//close
	destBf.Close();
}
