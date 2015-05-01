#if !defined(__PART_MODEL_H__)
#define __PART_MODEL_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include "mt_log.h"
namespace mt {
// roughly pInfo in PLL
class PartModelInfo {
    public:
    std::size_t GetNumRateCategories() const {
        return numRateCategories;
    }
    std::size_t GetNumStates() const {
        return numStates;
    }
    void SetExecuteMask(bool v) {
        executeMask = v;
    }
    bool GetExecuteMask() const {
        return executeMask;
    }
    bool CorrectForAscBias() const {
        return correctForAscBias;
    }
    bool GetUseRecom() const {
        return useRecom;
    }
    private:
    bool executeMask;
    std::size_t numRateCategories;
    std::size_t numStates;
    bool correctForAscBias;
    bool useRecom;
};

using PartModelInfoVec = std::vector<PartModelInfo>;


} // namespace mt
#endif