#if !defined(__LOG_H__)
#define __LOG_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>

template<typename T>
inline void _debug_vecl(const char * file, const char * func, int lineno, const char * l, const T &v) {
    std::cerr << file << " line " << lineno << " func " << func << " " << l << " len = " << v.size() << " vect = ";
    for (auto i : v) {
        std::cerr << i << ' ';
    }
    std::cerr << std::endl;
}
template<typename T>
inline void _debug_val(const char * file, const char * func, int lineno, const char * l, const T &v) {
    std::cerr << file << " line " << lineno << " func " << func << " " << l << " = " << v << std::endl;
}

inline void _debug_cla(const char * file,
                const char * func,
                int lineno,
                const char * v,
                const double *cla,
                unsigned nRateCats,
                unsigned nStates,
                unsigned numChars) {
    std::cerr << file << " line " << lineno << " func " << func << " " << v;
    std::cerr << " cla nr=" << nRateCats << " ns=" << nStates << " nc=" << numChars <<"\n";
    for (auto i = 0U; i < numChars; ++i) {
        for (auto ri = 0U; ri < nRateCats; ++ri) {
            for (auto fromState = 0U ; fromState < nStates; ++fromState) {
                std::cerr << cla[i*nRateCats*nStates + ri*nStates + fromState] << " ";
            }
        }
        std::cerr << "\n";
    }
    std::cerr << "/cla\n";
}


#define _DEBUG_VEC(a) (_debug_vecl(__FILE__, __FUNCTION__, __LINE__, "" #a "", a))
#define _DEBUG_VAL(a) (_debug_val(__FILE__, __FUNCTION__, __LINE__, "" #a "", a))
#define _DEBUG_CLA(a, r, s, c) (_debug_cla(__FILE__, __FUNCTION__, __LINE__, "" #a "", a, r, s, c))

#endif
