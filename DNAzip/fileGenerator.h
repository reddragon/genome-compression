#ifndef FILEGENERATOR_H_
#define FILEGENERATOR_H_
#include <string>
using namespace std;

enum opType { SNP, DELETION, INSERTION };

class fileGenerator
 { 
 public:     
     virtual ~fileGenerator() {}
     //This forces every derived class to have 
     //its own generate function
     virtual void generate()=0;      
 };

#endif /*FILEGENERATOR_H_*/
