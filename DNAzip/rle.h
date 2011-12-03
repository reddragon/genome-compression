#ifndef _RLE_H
#define _RLE_H

char digit2char(int digit);

string num2base62(int num);

string bitvec2compress(bool *vec, int len);
#endif
