#if !defined(__PART_MODEL_H__)
#define __PART_MODEL_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include "mt_log.h"

namespace mt {
class MTInstance;
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
    PartModelInfo(MTInstance & inst)
        :instance(inst) {
    }
    // delegation to instance for things, that (logically) could be part-specific, but currently aren't
    bool GetUseRecom() const;
    private:
    const MTInstance &instance;
    bool executeMask;
    std::size_t numRateCategories;
    std::size_t numStates;
    bool correctForAscBias;
};

using PartModelInfoVec = std::vector<PartModelInfo>;


} // namespace mt
#endif