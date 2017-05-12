#ifndef GENERAL_H
#define GENERAL_H

#include <stdlib.h>
#include <stdint.h>

#if !defined(max) && !defined(min)
size_t max(size_t a, size_t b);
size_t min(size_t a, size_t b);
#endif
int ispow2(size_t n);
int ispow2x(uintmax_t n);
int lg2(size_t n);
int lg2x(uintmax_t n);

#endif // GENERAL_H

