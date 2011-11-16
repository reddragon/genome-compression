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

#define MAXLEN 1100000
vector<bool> getBits(string chr)
{
	ifstream ifs;
	ifs.open(("bitvector/" + chr).c_str());
	char * bigline = (char *) malloc(MAXLEN), hbuf[200];
	assert(bigline);
	ifs.getline(hbuf, 200); /* Getting the Header */
	ifs.getline(bigline, MAXLEN);
	ifs.close();
	
	vector<bool> bits;
	int len = strlen(bigline);
	for(int i = 0; i < len; i++)
		bits.push_back(bigline[i]-'0');
	free(bigline);
	assert(len == bits.size());
	return bits;
}

vector<bool> getVecBool(int n, int min_len = 0)
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

void prependVec(vector<bool> & dest, vector<bool> src)
{
	dest.insert(dest.begin(), src.begin(), src.end());
}

void appendVec(vector<bool> & dest, vector<bool> src)
{
	dest.insert(dest.end(), src.begin(), src.end());
}



int toNum(vector<bool> vec)
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

void printVec(vector<bool> vec, bool new_line = 0)
{
	for(int i = 0; i < vec.size(); i++)
		cout << vec[i];
	if(new_line)
		cout << endl;
}

vector<bool> fetchVec(vector<bool> & dest, int len)
{
	vector<bool> ans = vector<bool>(dest.begin(), dest.begin() + len);
	dest.erase(dest.begin(), dest.begin() + len);
	return ans;
}

void decode(vector<bool> & encoded_str, vector<bool> & decoded_str, map< vector<bool>, vector<bool> > encoding)
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

void huffmanEncode(vector<bool> & bits, vector<bool> & encoded_str, int kmer_len = 5)
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

void huffmanDecode(vector<bool> &encoded_str, vector<bool> &decoded_str, int kmer_len = 5)
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

int osz, esz;
void doHuffman(vector<bool> bits, int kmer_len = 5)
{
	vector<bool> encoded_str, decoded_str;
	encoded_str.clear();
	huffmanEncode(bits, encoded_str, kmer_len);
	cout << "Compression Ratio: " << (encoded_str.size() * 1.0 / (bits.size())) << " ";
	huffmanDecode(encoded_str, decoded_str, kmer_len);
	assert(bits == decoded_str);
	osz += bits.size(); esz += encoded_str.size();
	cout << "Pass" << endl;
}

int main()
{
	/*
	vector<bool> test;
	test.push_back(1); test.push_back(0); test.push_back(1);
	test.push_back(1); test.push_back(0); test.push_back(1);
	test.push_back(1); test.push_back(0); test.push_back(1);
	cout << toNum(fetchVec(test,3)) << " ";
	cout << test.size() << endl;
	*/
	
	osz = 0; esz = 0;
	for(int i = 0; i <= 24; i++)
	{
		string chr = "chromosome";
		chr.push_back((i/10)+'0');
		chr.push_back((i%10)+'0');
		cout << chr << endl;
		doHuffman(getBits(chr));
	} 
	cout << "Original Size: " << osz << " " << "Compressed Size: " << esz << " " << "Compression: " << (esz*1.0)/osz << endl;
}
