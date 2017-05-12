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

uintmax_t large_rand()
{
    const size_t max_bits = lg2x(UINTMAX_MAX)+1, rand_bits = lg2(RAND_MAX)+1;
    size_t i = 0;
    uintmax_t value = 0;

    for (; i < (max_bits+rand_bits-1)/rand_bits; i += rand_bits)
    {
        value <<= rand_bits;
        value |= rand();
    }

    return value;
}

void memswap(void *a, void *b, size_t length)
{
    unsigned char *pa = a, *pb = b, c;
    while (length--)
    {
        c = *pa;
        *pa++ = *pb;
        *pb++ = c;
    }
}

void memreverse(void *dst, size_t elementsize, size_t length)
{
    size_t lo = 0, hi = length*elementsize;
    while (lo + elementsize < hi)
    {
        memswap(dst + lo, dst + hi - elementsize, elementsize);
        lo += elementsize;
        hi -= elementsize;
    }
}

#if defined(__linux__) || defined(__unix__)
intmax_t time_ms(intmax_t *tp)
{
    struct timeval tv;
    intmax_t result;

    gettimeofday(&tv, NULL);
    result = tv.tv_sec;
    result *= 1000;
    result += tv.tv_usec / 1000;
    if (tp != NULL)
        *tp = result;
    return result;
}
#else
intmax_t time_ms(intmax_t *tp)
{
    intmax_t result;

    result = time(NULL);
    result *= 1000;
    if (tp != NULL)
        *tp = result;
    return result;
}
#endif

#ifndef SLEEP_DEFINED
unsigned int sleep(unsigned int seconds)
{
    Sleep(seconds * 1000);
    return 0;
}
#define SLEEP_DEFINED
#endif
