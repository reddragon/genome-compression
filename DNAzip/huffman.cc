/* Copyright (c) 2008 the authors listed at the following URL, and/or
the authors of referenced articles or incorporated external code:
http://en.literateprograms.org/Huffman_coding_(C_Plus_Plus)?action=history&offset=20070405195734

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Huffman_coding_(C_Plus_Plus)?oldid=9682
*/

#include "huffman.h"
#include <map>

#include <iostream>
#include <ostream>

#include <algorithm>
#include <iterator>
#include <string>


std::ostream& operator<<(std::ostream& os, std::vector<bool> vec)
{
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<bool>(os, ""));
  return os;
}




/*int main()
{
  std::map<int, double> frequencies;
  frequencies[1] = 0.124167;
  frequencies[2] = 0.0820011;
  frequencies[3] = 0.0969225;
  frequencies[4] = 0.0768052;
  frequencies[5] = 0.0764055;
  frequencies[6] = 0.0714095;
  frequencies[7] = 0.0706768;
  frequencies[8] = 0.0668132;
  frequencies[9] = 0.0448308;
  frequencies[10] = 0.0363709;
  frequencies[11] = 0.0350386;
  frequencies[12] = 0.0344391;
  frequencies[13] = 0.028777;
  frequencies[14] = 0.0281775;
  frequencies[15] = 0.0235145;
  frequencies[16] = 0.0203171;
  frequencies[17] = 0.0189182;
  frequencies[18] = 0.0181188;
  frequencies[19] = 0.0135225;
  frequencies[20] = 0.0124567;
  frequencies[21] = 0.0106581;
  frequencies[22] = 0.00393019;
  frequencies[23] = 0.00219824;
  frequencies[24] = 0.0019984;
  frequencies[25] = 0.0009325;
  frequencies[26] = 0.000599;
  
  Hufftree<int, double> hufftree(frequencies.begin(), frequencies.end());

  for( int ch = 1; ch <= 29; ++ch )
  {
    std::cout << ch << ": " << hufftree.encode(ch) << "\n";
  }
  
  std::vector<bool> encoded = hufftree.encode(29);
  std::vector<int> destination;
  hufftree.decode(encoded, std::back_inserter(destination));
  
  std::vector<int>::iterator it;
  for ( it = destination.begin(); it!= destination.end(); ++it )
  	 cout << *it << " ";      
}*/

  /*std::string source = "literateprogramming";

  std::cout << source << "\n";

  std::vector<bool> encoded = hufftree.encode(source.begin(), source.end());
  std::cout << encoded << "\n";

  std::string destination;
  hufftree.decode(encoded, std::back_inserter(destination));
  std::cout << destination << "\n";*/
//}
//}
