#if !defined(__PART_DATA_H__)
#define __PART_DATA_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include "mt_log.h"
#include "mt_instance.h"

namespace mt {
class MTInstance;
class TraversalInfo;
class SlotIndices {
    public:
        const int p;
        const int q;
        const int r;
        SlotIndices(const MTInstance &instance, const TraversalInfo & tInfo);
};

class PartCalcs {
    public:
    double * GetCLArray(std::size_t slotNumber) {
        return &(internalCLVector.at(slotNumber).at(0));
    }
    private:
    std::vector<std::vector<double> > internalCLVector;
};

// part of  pInfo in PLL
class PartData {
    mutable PartCalcs partCalcs;
    public:
    std::size_t GetNumStates() const {
        return numStates;
    }
    double * GetCLArray(std::size_t slotNumber) const {
        return partCalcs.GetCLArray(slotNumber);
    }
    unsigned char * GetYArrayPtr(std::size_t yNumber) const;
    double * GetAscArrayPtr(std::size_t n) const; // pointer math... * pr->partitionData[model]->ascOffset
    unsigned int * GetGapArrayPtr(std::size_t n) const; // pointer math... * pr->partitionData[model]->gapVectorLength
    double * GetGapColumn(std::size_t n) const; // * states * rateHet];
    private:
    bool executeMask;
    std::size_t numRateCategories;
    std::size_t numStates;
    bool correctForAscBias;
    bool useRecom;
};
} // namespace mt
#endif