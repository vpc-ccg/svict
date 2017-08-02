
#ifndef Record_H
#define Record_H

#include <string>
#include <cstdlib>
#include "common.h"

const int MAXLEN = 8 * MB;

enum StringField {
    RN,
    CHR,
    CIGAR,
    P_CHR,
    SEQ,
    QUAL,
    OPT
};
enum IntField {
    MF,
	LOC,
    MQ,
    P_LOC,
    TLEN
};

class Record {
    char line[MAXLEN];
    char *strFields[7];
    int32_t intFields[5];

private:
    friend class BAMParser;
    friend class SAMParser;

public:
    const char* getReadName() const { return strFields[RN]; }
    const int getMappingFlag() const { return intFields[MF]; }
    const char* getChromosome() const { return strFields[CHR]; }
    const size_t getLocation() const { return intFields[LOC] - 1; }
    const char getMappingQuality() const { return intFields[MQ]; }
    const char* getCigar() const { return strFields[CIGAR]; }
    const char* getPairChromosome() const { return strFields[P_CHR]; }
    const size_t getPairLocation() const { return intFields[P_LOC] - 1; }
    const int getTemplateLength() const { return intFields[TLEN]; }
    const char* getSequence() const { return strFields[SEQ]; }
    const char* getQuality() const { return strFields[QUAL]; }
    const char* getOptional() const { return strFields[OPT]; }

    std::string getFullRecord() const {
        return S(
            "%s %d %s %lu %d %s %s %lu %d %s %s p",
            getReadName(),
            getMappingFlag(),
            getChromosome(),
            getLocation(),
            getMappingQuality(),
            getCigar(),
            getPairChromosome(),
            getPairLocation(),
            getTemplateLength(),
            getSequence(),
            getQuality(),
            getOptional()
        );
    }

    void testRecords() const {
        LOG(
            "%s %d %s %lu %d %s %s %lu %d %s %s %s\n",
            getReadName(),
            getMappingFlag(),
            getChromosome(),
            getLocation(),
            getMappingQuality(),
            getCigar(),
            getPairChromosome(),
            getPairLocation(),
            getTemplateLength(),
            getSequence(),
            getQuality(),
            getOptional()
        );
    }
};

#endif // RECORD_H
