#ifndef BITHUFFMAN_H_INC
#define BITHUFFMAN_H_INC

#include <iostream>
#include "huffman.h"
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>

using namespace std;

class BitHuffman {
    
    public:

    static vector<bool> getVecBool(int n, int min_len);
    static void prependVec(vector<bool> & dest, vector<bool> src);
    static void appendVec(vector<bool> & dest, vector<bool> src);
    static int toNum(vector<bool> vec);
    static void printVec(vector<bool> vec, bool new_line);
    static vector<bool> fetchVec(vector<bool> & dest, int len);
    static void decode(vector<bool> & encoded_str, \
		vector<bool> & decoded_str, map< vector<bool>, \
		vector<bool> > encoding);
	static void huffmanEncode(vector<bool> & bits, \
		vector<bool> & encoded_str, int kmer_len);
	static void huffmanDecode(vector<bool> &encoded_str, \
		vector<bool> &decoded_str, int kmer_len);
};

vector<bool> BitHuffman::getVecBool(int n, int min_len = 0)
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

void BitHuffman::prependVec(vector<bool> & dest, vector<bool> src)
{
	dest.insert(dest.begin(), src.begin(), src.end());
}

void BitHuffman::appendVec(vector<bool> & dest, vector<bool> src)
{
	dest.insert(dest.end(), src.begin(), src.end());
}

int BitHuffman::toNum(vector<bool> vec)
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

void BitHuffman::printVec(vector<bool> vec, bool new_line = 0)
{
	for(int i = 0; i < vec.size(); i++)
		cout << vec[i];
	if(new_line)
		cout << endl;
}

vector<bool> BitHuffman::fetchVec(vector<bool> & dest, int len)
{
	vector<bool> ans = vector<bool>(dest.begin(), dest.begin() + len);
	dest.erase(dest.begin(), dest.begin() + len);
	return ans;
}

void BitHuffman::decode(vector<bool> & encoded_str, vector<bool> & decoded_str, map< vector<bool>, vector<bool> > encoding)
{
	vector<bool> temp;

	int c = 0;

	for(int i = 0; i < encoded_str.size(); i++)
	{
		temp.push_back(encoded_str[i]);
		if(encoding.find(temp) != encoding.end())
		{
			appendVec(decoded_str, encoding[temp]);
			temp.clear();
			continue;
		}
	}
	if(temp.size() != 0)
		appendVec(decoded_str, encoding[temp]);
}

void BitHuffman::huffmanEncode(vector<bool> & bits, vector<bool> & encoded_str, int kmer_len = 5)
{
	vector< vector<bool> > kmer_strs;
	vector<bool> kmer;
	map< vector<bool>, int> counts;
	
	for(int i = 0; i < bits.size(); i++)
	{
		if(i && (i % kmer_len == 0)) 
		{ 
			kmer_strs.push_back(kmer); 
			counts[kmer] += 1; 
			kmer.clear(); 
		} 
		kmer.push_back(bits[i]);
	}

	if(kmer.size() != 0) 
	{ 
		while(kmer.size() != kmer_len) 
			kmer.push_back(0);
		kmer_strs.push_back(kmer);
		counts[kmer] += 1; 
		kmer.clear();
	}
	Hufftree<vector<bool>, int> hufftree(counts.begin(), counts.end());
	
	encoded_str = hufftree.encode(kmer_strs.begin(), kmer_strs.end());
	vector<bool> dict;

	for(int i = 0; i < (1<<kmer_len); i++)
	{
		vector<bool> vec = getVecBool(i, kmer_len);
		if(counts.find(vec) == counts.end())
		{
			appendVec(dict, getVecBool(0,4));			
		}
		else
		{
			vector<bool> enc_vec = hufftree.encode(getVecBool(i,kmer_len));
			appendVec(dict, getVecBool(enc_vec.size(), 4));
			appendVec(dict, enc_vec);
		}
	}
	
	// Sanity check for the dict	
	vector< bool > c_dict = dict;
	
	for(int i = 0; i < (1<<kmer_len); i++)
	{
		vector<bool> vec = getVecBool(i,kmer_len);
		int len = toNum(fetchVec(c_dict, 4));
		
			
		if(len == 0)
		{
			assert(counts.find(vec) == counts.end());
			continue;
		}
		vector<bool> rep = fetchVec(c_dict, len);
		//encoding[rep] = vec;
		assert(hufftree.encode(vec) == rep); 
	}  
	
	//cout << "dict is sane" << endl;
	// Prepending the number of padding we have done
	prependVec(encoded_str, getVecBool((kmer_strs.size() * kmer_len) - bits.size(), 3));
	//cout << "Padding: " << (kmer_strs.size() * kmer_len) - bits.size() << endl;
	vector<bool> dict_sz = getVecBool(dict.size(), 10);
	prependVec(encoded_str, dict);
	prependVec(encoded_str, dict_sz);
}

void BitHuffman::huffmanDecode(vector<bool> &encoded_str, vector<bool> &decoded_str, int kmer_len = 5)
{
	//cout << "Decompression begins" << endl;
	int dict_size = toNum(fetchVec(encoded_str, 10));
	//cout << "Dict Size: " << dict_size << endl;
	vector<bool> new_dict = fetchVec(encoded_str, dict_size);
	
	map< vector<bool>, vector<bool> > encoding;
	for(int i = 0; i < (1<<kmer_len); i++)
	{
		vector<bool> vec = getVecBool(i,kmer_len);
		int len = toNum(fetchVec(new_dict, 4));
			
		vector<bool> rep = fetchVec(new_dict, len);
		encoding[rep] = vec;
	}  

	int bits_padding = toNum(fetchVec(encoded_str, 3));
	//cout << "Padding: " << bits_padding << endl;

	
	decode(encoded_str, decoded_str, encoding);
	//printVec(decoded_string,1);
	decoded_str.erase(decoded_str.end() - bits_padding, decoded_str.end());
	
	//cout << decoded_str.size() << endl;
}

#endif
