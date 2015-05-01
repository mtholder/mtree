#if !defined(__PART_DATA_H__)
#define __PART_DATA_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include "mt_log.h"
namespace mt {

class SlotIndices {
    public:
        const int p;
        const int q;
        const int r;
        SlotIndices(const MTInstance &instance, const TraversalInfo & tInfo) 
            :p(instance.GetUseRecom() ? tInfo.slot_p : tInfo.pNumber - instance.GetMxTips() - 1),
            q(instance.GetUseRecom() ? tInfo.slot_p : tInfo.pNumber - instance.GetMxTips() - 1),
            r(instance.GetUseRecom() ? tInfo.slot_p : tInfo.pNumber - instance.GetMxTips() - 1) {
        }
};
class PartCalcs {
    double * GetCLArray(std::size_t slotNumber) {
        &(internalCLVector.at(slotNumber)[0]);
    }
    private:
    std::vector<std::vector<double> > internalCLVector;
};

// part of  pInfo in PLL
class PartData {
    mutable PartCalcs partCalcs;
    public:
    allocCLArray
    std::size_t GetNumStates() const {
        return numStates;
    }
    double * const GetCLArray(std::size_t slotNumber) const {
        return partCalcs.GetCLArray(slotNumber)
v;    }
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