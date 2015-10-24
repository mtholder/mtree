#if !defined(__PARSIMONY__)
#define __PARSIMONY__


#include "mt_log.h"
#include "mt_instance.h"


namespace mt {

typedef unsigned char BitField;
typedef std::vector<BitField> BitFieldRow;

class ParsInfo {
    public:
        ParsInfo()
            :downPass(0),
            allSeen(0),
            score(0),
            numPatterns(0) {
            }
        std::size_t size() const {
            return numPatterns;
        }
        void calculateForTip(const BitFieldRow & data, MTInstance & instance);
        void calculateForInternal(ParsInfo & leftData, ParsInfo & rightData);
        void write(std::ostream & o) const;
        const BitField * downPass; // alias
        const BitField * allSeen; // alias
        const std::size_t * score; // alias

    private:
        std::size_t numPatterns;
        BitFieldRow downPassOwned;
        BitFieldRow allSeenOwned;
        std::vector<std::size_t> scoreOwned;

};

} //namespace

#endif
