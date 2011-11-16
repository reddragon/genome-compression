#include "input.h"

//read a vector of bits
void readBits( bit_file_c& bf, unsigned count, vector<bool>* bits )
{
	  /* read single bits */
	  int value = 0;
	  for ( unsigned i = 0; i < count; i++ )
	  {
	      value = bf.GetBit();
	      if(value == EOF)
	      {
	         cerr << "Error: reading bit" << endl;
	         bf.Close();
	         return;
	      }
	      if(value == 1)
	      {
	    	  bits->push_back( true );	    	  
	      }
	      else
	      {
	    	  bits->push_back( false);	    	  
	      }	      	    	  
	   }	
}

//read an array of bits
void readBitArrays( bit_file_c& bf, bool* bits, unsigned cnt )
{
	/* read single bits */
	int value = 0;
	for ( unsigned i = 0; i < cnt; i++ )
	{
	    value = bf.GetBit();
	    if(value == EOF)
	    {
		    cerr << "Error: reading bit" << endl;
		    bf.Close();
		    return;
		}
		if(value == 1)
		{
		   bits[i] = true;		   
		}
		else
		{
		   bits[i] = false;		   
		}
     }	
}

//read a vector of gene letters
void readBitGens( bit_file_c& bf, unsigned count, vector<char>* gens )
{				
	int value1 = 0;
	int value2 = 0;
	for( unsigned i = 0; i < count; i++ )
	{
	   value1 = bf.GetBit();
	   value2 = bf.GetBit();
				
	   if( value1 == 0 && value2 == 0 )
	   {
		  gens->push_back('A');
	   }
	   else if( value1 == 0 && value2 == 1 )
	   {
		  gens->push_back('T');
	   }
	   else if( value1 == 1 && value2 == 0 )
	   {
		  gens->push_back('C');
	   }
	   else
	   {
		  gens->push_back('G');
	   }						
	}	
}

//read a VINT 
unsigned readBitVINT( bit_file_c& bf )
{
	char b = (char)bf.GetBit();
	for( unsigned j = 1; j < 8; j++ )
	{
		b <<= 1;
		int value = bf.GetBit();
		b |= value;
	}
				
	signed i = b & 0x7F;	
	for( signed shift = 7; (b & 0x80) != 0; shift += 7 ) {
		  b = (char)bf.GetBit();
		  for( unsigned j = 1; j < 8; j++ )
		  {
			  b <<= 1;
			  int value = bf.GetBit();
			  b |= value;
		  }
	      i |= (b & 0x7F) << shift;
	}
	return (unsigned)i;	
}

//read a string
string readString( bit_file_c& bf )
{
	 
	int strSZ = bf.GetChar();	
	
	if( strSZ == EOF )
		return "";		
	char* str = new char[strSZ+1];
	for( int i = 0; i < strSZ; i++ )
	{
		str[i] = (char)bf.GetChar();
	}
	str[strSZ] = '\0';		
	return string(str);
}
