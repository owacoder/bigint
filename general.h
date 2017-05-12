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
uintmax_t large_rand();
void memswap(void *a, void *b, size_t length);
void memreverse(void *dst, size_t elementsize, size_t length);

/* TODO: time_ms() function for Windows */
#if defined(__linux__) || defined(__unix__)
# include <unistd.h>
# include <time.h>
# include <sys/time.h>
# define SLEEP_DEFINED
#else
# include <windows.h>
unsigned int sleep(unsigned int seconds);
#endif

intmax_t time_ms(intmax_t *tp);

#endif // GENERAL_H

