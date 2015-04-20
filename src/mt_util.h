#if !defined(__UTIL_H__)
#define __UTIL_H__

namespace mt {

template <typename T>
inline void enforce_bound(T & val, const T & minBound, const T & maxBound) {
    if (val < minBound) {
        val = minBound;
    } else if (val > maxBound) {
        val = maxBound;
    }
}

template <typename T>
inline void enforce_min_bound(T & val, const T & minBound) {
    if (val < minBound) {
        val = minBound;
    }
}

template <typename T>
inline void enforce_max_bound(T & val, const T & maxBound) {
    if (val > maxBound) {
        val = maxBound;
    }
}
template <typename T>
inline bool all_true(const T & vContainer) {
    for (const auto & v : vContainer) {
        if (!v) {
            return false;
        }
    }
    return true;
}



} //namespace mt

#endif
