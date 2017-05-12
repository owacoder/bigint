#ifndef GENERAL_H
#define GENERAL_H

#include <stdlib.h>
#include <stdint.h>

size_t max(size_t a, size_t b);
size_t min(size_t a, size_t b);
int ispow2(size_t n);
int ispow2x(uintmax_t n);
int lg2(size_t n);
int lg2x(uintmax_t n);

#endif // GENERAL_H

