#include "general.h"
#include <limits.h>

#if !defined(max) && !defined(min)
size_t max(size_t a, size_t b)
{
    return a > b? a: b;
}

size_t min(size_t a, size_t b)
{
    return a < b? a: b;
}
#endif

int ispow2(size_t n)
{
    if (!n) return 0;
    return (n & (n-1)) == 0;
}

int ispow2x(uintmax_t n)
{
    if (!n) return 0;
    return (n & (n-1)) == 0;
}

int lg2(size_t n)
{
    int l = -1;
    for (; n; n >>= 1)
        ++l;
    return l;
}

int lg2x(uintmax_t n)
{
    int l = -1;
    for (; n; n >>= 1)
        ++l;
    return l;
}
