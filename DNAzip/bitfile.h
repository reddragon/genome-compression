#ifndef BITFILE_H_
#define BITFILE_H_
#include <iostream>
#include <fstream>

/***************************************************************************
*                            TYPE DEFINITIONS
***************************************************************************/
typedef enum
{
    BF_READ = 0,
    BF_WRITE = 1,
    BF_APPEND= 2,
    BF_NO_MODE
} BF_MODES;

typedef enum
{
    BF_UNKNOWN_ENDIAN,
    BF_LITTLE_ENDIAN,
    BF_BIG_ENDIAN
} endian_t;

class bit_file_c
{
    public:
        bit_file_c(void);
        bit_file_c(const char *fileName, const BF_MODES mode);
        virtual ~bit_file_c(void);

        /* open/close bit file */
        void Open(const char *fileName, const BF_MODES mode);
        void Close(void);

        /* toss spare bits and byte align file */
        int ByteAlign(void);

        /* fill byte with ones or zeros and write out results */
        int FlushOutput(const unsigned char onesFill);

        /* get/put character */
        int GetChar(void);
        int PutChar(const int c);

        /* get/put single bit */
        int GetBit(void);
        int PutBit(const int c);

        /* get/put number of bits */
        int GetBits(void *bits, const unsigned int count);
        int PutBits(void *bits, const unsigned int count);

        /* get/put number of bits to/from integer types (short, int, ...)*/
        /* size is size of data structure pointed to by bits.            */
        int GetBitsInt(void *bits, const unsigned int count,
            const size_t size);
        int PutBitsInt(void *bits, const unsigned int count,
            const size_t size);

        /* status */
        bool eof(void);
        bool good(void);
        bool bad(void);

    private:
        std::ifstream *m_InStream;      /* input file stream pointer */
        std::ofstream *m_OutStream;     /* output file stream pointer */
        endian_t m_endian;              /* endianess of architecture */
        char m_BitBuffer;               /* bits waiting to be read/written */
        unsigned char m_BitCount;       /* number of bits in bitBuffer */
        BF_MODES m_Mode;                /* open for read, write, or append */

        /* endianess aware methods used by GetBitsInt/PutBitsInt */
        int GetBitsLE(void *bits, const unsigned int count);
        int PutBitsLE(void *bits, const unsigned int count);

        int GetBitsBE(void *bits, const unsigned int count,
            const size_t size);
        int PutBitsBE(void *bits, const unsigned int count,
            const size_t size);
};

#endif  /* ndef BITFILE_H_ */
