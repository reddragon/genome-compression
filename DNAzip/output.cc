#include "output.h"

//write a vector of bits
void writeBits( bit_file_c& bf, vector<bool>* bits )
{
	for (unsigned i = 0; i < bits->size(); i++)
	{	 
		 signed value = 0;
		 if( (*bits)[i] )
			value = 1; 
	     if(bf.PutBit(value) == EOF)
	     {
	          cerr << "Error: writing bit" << endl;
	          bf.Close();
	          return ;
	     }	        
	}
}

//write an array of bits 
void writeBitArrays( bit_file_c& bf, bool* bits, unsigned cnt )
{
	for (unsigned i = 0; i < cnt; i++)
	{	 
		signed value = 0;			
		
		if( bits[i] )		
			value = 1;			
		
		if(bf.PutBit(value) == EOF)
		{
		    cerr << "Error: writing bit" << endl;
		    bf.Close();
		    return ;
		}	        
	}
}

//write a vector of gene letters
void writeBitGens( bit_file_c& bf, vector<char> * gens )
{
	vector<char>::iterator iter = gens->begin();		
	while( iter != gens->end() ) {
	char gen = *iter;	    
	    
	if( gen == 'a' || gen == 'A' ) //00: A
	{ 
		bf.PutBit(0);
	    bf.PutBit(0);
	}
	else if( gen == 't' || gen == 'T' ) //01: T
	{
		bf.PutBit(0);
		bf.PutBit(1);
	}
	else if( gen == 'c' || gen == 'C' ) //10: C
	{
		bf.PutBit(1);
		bf.PutBit(0);
	}
	else  //11: G
	{
		bf.PutBit(1);
		bf.PutBit(1);
	}	    	    
	iter++;	
    }	
}

//write a VINT 
void writeBitVINT( bit_file_c& bf, unsigned num)
{
	int word_length = 5;
	vector<bool> * bits = new vector<bool> ();
	char mask = (1 << (word_length - 1)) - 1;
	do
	{
		if(num > mask) bits->push_back(1);
		else bits->push_back(0);
		
		char word = num & mask;
		//cout << (int)word << endl;
		
		for(int i = (word_length - 2); i >= 0; i--)
			bits->push_back((word & (1<<i)));
		
		num >>= (word_length - 1);
	} while(num);
	
	
	writeBits( bf, bits );
	//cout << num << " " << c << endl;
	delete bits;	
}


/* 
void writeBitVINT( bit_file_c& bf, unsigned num )
{
	vector<bool>* bits = new vector<bool>();
	int mask, c = 0;
	unsigned x = num;
	char chrByte;
	while ( num >= 128 ) {
		 chrByte = (char)((num & 0x7f) | 0x80);		     
		 mask = 128;
		 for( unsigned i = 0; i < 8; i++ )
		 {
		 	 if( chrByte & mask )
		    	 bits->push_back( true );
		     else
		    	 bits->push_back( false );
		     mask >>= 1;
		  }		     
		  num >>= 7;
		  //c++;
	}
	
	chrByte = (char)num;
	mask = 128;
	for( unsigned i = 0; i < 8; i++ )
	{
		if( chrByte & mask )
		   bits->push_back( true );
		else
		   bits->push_back( false );
		mask >>= 1;
		
	}
	
	writeBits( bf, bits );
	cerr << x << endl;
	//cout << num << " " << c << endl;
	delete bits;	
}*/

//write a string
void writeString( bit_file_c& bf, string& str )
{
	int strSZ = str.size();	
	bf.PutChar( strSZ );
	for( int i = 0; i < strSZ; i++ )
	{
		bf.PutChar( str[i] );
	}
}

