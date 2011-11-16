#include "watsonFileGenerator.h"

watsonFileGenerator::watsonFileGenerator(const string& SNPFile,
        const string& indelFile, 
        const string& destFile):SNPFile(SNPFile),indelFile(indelFile),destFile(destFile)
{
}

void watsonFileGenerator::generate()
{
	string chromosome;
	string pos;
	string ref_alt;
	
	//destion file with unified format
	ofstream destStream( destFile.c_str() );
	//convert SNP file into destination file
	ifstream snpStream( SNPFile.c_str() );
	string line;
	getline( snpStream, line);	
	string prevChromosome;
	string prevPos;
	
	while( getline( snpStream, line ) )
	{		
		char* token;		
		token = strtok ( (char *)line.c_str(), "," ); //ID
		token = strtok ( NULL, "," ); //Chromosome
		chromosome = token;		
		token = strtok ( NULL, ","); //Start
		pos = token;		
		token = strtok ( NULL, "," ); //Ref/Alt
		ref_alt = token;
		
		if( chromosome.compare(prevChromosome) == 0 && pos.compare(prevPos) == 0 )
			continue;		
		destStream<<SNP<<","<<chromosome<<","<<pos<<","<<ref_alt<<endl;
		
		prevChromosome = chromosome;
		prevPos = pos;
	}
	snpStream.close();
	//convert insertion/deletion file into destination file
	ifstream indelStream( indelFile.c_str() );	
	getline( indelStream, line);			
	while( getline( indelStream, line ) )
	{		
		char* token;		
		token = strtok( (char *)line.c_str(), "," ); //ID
		token = strtok( NULL, "," ); //Chromosome
		chromosome = token;		
		token = strtok ( NULL, ","); //Start
		pos = token;
		token = strtok ( NULL, ","); //End
		token = strtok ( NULL, "," ); //Ref/Alt
		ref_alt = token;	
		
		if( chromosome.compare(prevChromosome) == 0 && pos.compare(prevPos) == 0 )
					continue;	
		
		if( ref_alt.find("-/") != string::npos )
			destStream<<INSERTION<<","<<chromosome<<","<<pos<<","<<ref_alt<<endl;
		else
			destStream<<DELETION<<","<<chromosome<<","<<pos<<","<<ref_alt<<endl;
		
		prevChromosome = chromosome;
		prevPos = pos;		
	}
	indelStream.close();
	destStream.close();		
}
