# $Id$	
# Date: 09/26/2007
# Author: Yiming Lu

CC = g++
#CPPFLAGS = -g -ggdb  -O0
CPPFLAGS = -O3

all: perftest

input.o: input.h 

output.o: output.h

watsonFileGenerator.o: watsonFileGenerator.h

kMerFreqGenerator.o: kMerFreqGenerator.h

dbSNPCompression.o: dbSNPCompression.h 

dbSNPDeCompression.o: dbSNPDeCompression.h

huffman.o: huffman.h

bitfile.o: bitfile.h

perftest.o: input.h output.h dbSNPCompression.h dbSNPDeCompression.h watsonFileGenerator.h kMerFreqGenerator.h 

perftest: perftest.o bitfile.o huffman.o input.o output.o dbSNPCompression.o dbSNPDeCompression.o kMerFreqGenerator.o watsonFileGenerator.o

clean:
	@rm perftest *.o
