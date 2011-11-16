#include "kMerFreqGenerator.h"

kMerFreqGenerator::kMerFreqGenerator(
		int kMer,
		const string& srcFile,
        const string& destFile
        ):kMer(kMer), srcFile(srcFile), destFile(destFile)
{				
}

void kMerFreqGenerator::generate()
{
	ifstream srcStream( srcFile.c_str() );
	string line;
	getline( srcStream, line);			
				
	bit_file_c bf; 
	bf.Open( destFile.c_str(), BF_WRITE );		
								
	std::map<string, int> strFreq; //frequencies
					
	while( getline( srcStream, line ) )
	{						
		//Tokenize the line						
		char* token;
		token = strtok ( (char *)line.c_str(), "," ); //type
		token = strtok ( NULL, "," ); //chromsomeID
		token = strtok ( NULL, "," ); //start position				
		token = strtok ( NULL, "," ); //genLetters	
		string strToken = token;
		if( strToken.find("-/") == string::npos ) // if it is not insertion
			continue;
			
		int startPos = strToken.length()/2;		
		string genToken =  strToken.substr( startPos+1, startPos );
			
		signed numOfLongs = genToken.length()/kMer;
		for( signed i = 0; i < numOfLongs; i++ )
		{
			string subGenToken = genToken.substr(i*kMer, kMer);
					
			map<string,signed>::iterator iter;   				
			signed freq = 0;
					
			if( (iter = strFreq.find(subGenToken)) != strFreq.end() )
			{
				freq = (signed)iter->second;
			}
					
			freq++;
					
			strFreq[subGenToken] = freq;
		 }		
		  
		if( genToken.length()%kMer > 0 )
		{
			string subGenToken = genToken.substr( numOfLongs*kMer, genToken.length()%kMer );
								
			map<string,signed>::iterator iter;   				
			signed freq = 0;
								
			if( (iter = strFreq.find(subGenToken)) != strFreq.end() )
			{
				freq = (signed)iter->second;
			}
								
			freq++;
								
			strFreq[subGenToken] = freq;			
		}		
	   }		
		
		writeBitVINT(bf, strFreq.size());
		std::map<string, int>::iterator iter = strFreq.begin();
		while( iter != strFreq.end() )
		{
			//destionationStream <<  iter->first <<"," << iter->second <<endl;
			string gens = iter->first;
			writeString( bf,  gens );
			writeBitVINT(bf, iter->second);
			iter++;		
		}	
		
		bf.Close();
		srcStream.close();								
}
