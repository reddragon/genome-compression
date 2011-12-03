#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
using namespace std;
#include "rle.h"

char digit2char(int digit)
{
        if(digit >=62)
                exit(1);
        if(digit < 10) {
                return '0'+digit;
        }
        if(digit < 36) {
                return digit -10 + 'A';
        }
        return digit - 36 + 'a';
}

string num2base62(int num)
{
        string s;
        while(num > 0) {
                s.insert(s.end(), digit2char(num % 62));
                num = num / 62;
        }
        return string(s.rbegin(), s.rend());
}

string bitvec2compress(bool *vec, int len)
{
        string final;
        int count = 1;
        bool prev = vec[0];
        final.insert(final.end(), prev?'+':'-');
        for(int i=1;i<=len;i++) {
                bool next = vec[i];
                if(next == prev && i != len) {
                        count ++;
                        prev = next;
                        continue;
                }
                else {
                        //cout << count << " " << prev <<"s";
                        string ins = num2base62(count);
                        //cout << " "<< ins << endl;
                        int len = ins.length();
                        char precede = ' ';
                        switch(len) {
                                case 1: precede = ' '; break;
                                case 2: precede = '!'; break;
                                case 3: precede = '@'; break;
                                case 4: precede = '#'; break;
                                case 5: precede = '$'; break;
                                case 6: precede = '%'; break;
                                default: cout << "no";
                        }
                        if(len == 1)
                                final = final + ins;
                        else {
                                final.insert(final.end(), precede);
                                final = final + ins;
                        }
                        prev = next;
                        count = 1;
                }
        }
        return final;
}
