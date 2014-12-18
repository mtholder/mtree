#if !defined(__LOG_H__)
#define __LOG_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
namespace mt {

template<typename T>
void _debug_vec(const T &v) {
    std::cerr << "_debug_vec len = " << v.size() << " vect = ";
    for (auto i : v) {
        std::cerr << i << ' ';
    }
    std::cerr << std::endl;
}
template<typename T>
void _debug_vecl(const char * file, const char * func, int lineno, const char * l, const T &v) {
    std::cerr << file << " line " << lineno << " func " << func << " " << l << " len = " << v.size() << " vect = ";
    for (auto i : v) {
        std::cerr << i << ' ';
    }
    std::cerr << std::endl;
}
template<typename T>
void _debug_val(const char * l, const T &v) {
    std::cerr << "_debug_val " << l << " = " << v << std::endl;
}

#define _DEBUG_VEC(a) (_debug_vecl(__FILE__, __FUNCTION__, __LINE__, "" #a "", a))


} //namespace
#endif
