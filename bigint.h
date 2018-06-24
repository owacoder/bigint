#ifndef BIGINT_H
#define BIGINT_H

#include <memory.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>
#include "general.h"

#ifndef BIGINT_WORD_SIZE
#if SIZE_MAX == 0xffffffffffffffff // Assume it's a 64-bit platform
#define BIGINT_WORD_SIZE 32
#elif SIZE_MAX == 0xffffffff // Assume a 32-bit platform
#define BIGINT_WORD_SIZE 16
#elif SIZE_MAX == 0xffff // Assume a 16-bit platform
#define BIGINT_WORD_SIZE 8
#else
#error BIGINT_WORD_SIZE must be specified. Cannot detect word size of target platform.
#endif
#endif

#if BIGINT_WORD_SIZE == 64 && defined(__GNUC__)
// Size of leaves, or digit cells. Must be unsigned.
typedef uint64_t bi_leaf;
// Size of signed leafs. Must be signed. An unsigned leaf must be able to store all the positive values.
typedef int64_t bi_signed_leaf;
// Size of double-leaves, must fit the maximum leaf value squared. Must be unsigned.
typedef __uint128_t bi_dleaf;
// Size of largest unsigned integer values
typedef __uint128_t bi_uintmax;
// Size of largest signed integer values
typedef __int128_t bi_intmax;
#define BIGINT_LEAFSIZE (8)
#elif BIGINT_WORD_SIZE == 32
// Size of leaves, or digit cells. Must be unsigned.
typedef uint32_t bi_leaf;
// Size of signed leafs. Must be signed. An unsigned leaf must be able to store all the positive values.
typedef int32_t bi_signed_leaf;
// Size of double-leaves, must fit the maximum leaf value squared. Must be unsigned.
typedef uint64_t bi_dleaf;
// Size of unsigned maximum values
typedef uintmax_t bi_uintmax;
// Size of signed maximum values
typedef intmax_t bi_intmax;
#define BIGINT_LEAFSIZE (4)
#elif BIGINT_WORD_SIZE == 16
// Size of leaves, or digit cells. Must be unsigned.
typedef uint16_t bi_leaf;
// Size of signed leafs. Must be signed. An unsigned leaf must be able to store all the positive values.
typedef int16_t bi_signed_leaf;
// Size of double-leaves, must fit the maximum leaf value squared. Must be unsigned.
typedef uint32_t bi_dleaf;
// Size of unsigned maximum values
typedef uintmax_t bi_uintmax;
// Size of signed maximum values
typedef intmax_t bi_intmax;
#define BIGINT_LEAFSIZE (2)
#elif BIGINT_WORD_SIZE == 8
// Size of leaves, or digit cells. Must be unsigned.
typedef uint8_t bi_leaf;
// Size of signed leafs. Must be signed. An unsigned leaf must be able to store all the positive values.
typedef int8_t bi_signed_leaf;
// Size of double-leaves, must fit the maximum leaf value squared. Must be unsigned.
typedef uint16_t bi_dleaf;
// Size of unsigned maximum values
typedef uintmax_t bi_uintmax;
// Size of signed maximum values
typedef intmax_t bi_intmax;
#define BIGINT_LEAFSIZE (1)
#else
#error BIGINT_WORD_SIZE must be defined to one of the following: 8, 16, or 32 (or 64 on GCC targets)
#endif

// Minimum number of leaves allocated for a bigint. Must be greater than 0.
#ifndef BIGINT_MINLEAFS
#define BIGINT_MINLEAFS (8)
#endif

#define BIGINT_LEAFMAX ((bi_leaf) ~0ULL)
#define BIGINT_SIGNED_LEAFMAX ((bi_signed_leaf) (BIGINT_LEAFMAX >> 1))
#define BIGINT_LEAFBYTES (BIGINT_LEAFSIZE*CHAR_BIT/8)
#define BIGINT_LEAFBITS (BIGINT_LEAFSIZE*CHAR_BIT)
#define BIGINT_LEAFS_PER_BI_INTMAX (sizeof(bi_intmax)/BIGINT_LEAFSIZE)

#ifndef BIGINT_KARATSUBA_MIN_LEAFS
#define BIGINT_KARATSUBA_MIN_LEAFS (640/BIGINT_LEAFBITS)
#endif

#ifndef BIGINT_TOOM_COOK_MIN_LEAFS
#define BIGINT_TOOM_COOK_MIN_LEAFS (10000/BIGINT_LEAFBITS)
#endif

#if defined(__GNUC__) && !defined(BIGINT_DISABLE_PTHREADS)
#define BIGINT_ENABLE_PTHREADS
#endif

#if (defined(_WIN32) || defined(_WIN64)) && !defined(BIGINT_DISABLE_WINTHREADS)
#define BIGINT_ENABLE_WINTHREADS
#endif

#if defined(BIGINT_ENABLE_WINTHREADS)
#undef BIGINT_ENABLE_PTHREADS
#endif

#if !defined(BIGINT_DISABLE_LIBMATH)
#define BIGINT_ENABLE_LIBMATH
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    char *string;
    size_t len;
} bigint_string;

bigint_string bis_new(char *s, size_t len);
void bis_destroy(bigint_string s);

#define BIGINT_FLAG_SIGN     (unsigned long) 0x01
#define BIGINT_FLAG_SIGN_SHIFT 0
#define BIGINT_FLAG_GROWABLE (unsigned long) 0x02
#define BIGINT_FLAG_GROWABLE_SHIFT 1
#define BIGINT_FLAG_FREEABLE (unsigned long) 0x04
#define BIGINT_FLAG_FREEABLE_SHIFT 2
#define BIGINT_FLAG_DESTROYABLE (unsigned long) 0x08
#define BIGINT_FLAG_DESTROYABLE_SHIFT 3

typedef struct
{
    bi_leaf *data;
    size_t size;
    // Bit 0 is sign (unset is "positive")
    // Bit 1 is "growable" flag (unset is "not growable")
    // Bit 2 is "freeable" flag (unset is "data storage not freeable")
    // Bit 3 is "destroyable" flag (unset is "struct not destroyable")
    size_t flags;
} bigint;

bigint bi_zero();
bigint bi_one();
bigint bi_two();
bigint bi_minus_one();
bigint *bi_new();
bigint *bi_init(bigint *bi);
bigint *bi_new_valueu(bi_leaf value);
bigint *bi_new_value(bi_signed_leaf value);
bigint *bi_new_valuelu(bi_uintmax value);
bigint *bi_new_valuel(bi_intmax value);
bigint *bi_init_valueu(bigint *bi, bi_leaf value);
bigint *bi_init_value(bigint *bi, bi_signed_leaf value);
bigint *bi_init_valuelu(bigint *bi, bi_uintmax value);
bigint *bi_init_valuel(bigint *bi, bi_intmax value);
bigint *bi_copy(const bigint *bi);
bigint *bi_init_copy(bigint *dst, const bigint *src);
bigint *bi_copy_mag_to(bigint *dst, const bigint *src);
bigint *bi_copy_to(bigint *dst, const bigint *src);
bigint *bi_assignu(bigint *bi, bi_leaf value);
bigint *bi_assign(bigint *bi, bi_signed_leaf value);
bigint *bi_assignl(bigint *bi, bi_intmax value);
bigint *bi_assignlu(bigint *bi, bi_uintmax value);
bi_signed_leaf bi_bit(const bigint *bi, size_t bit);
int bi_set_bit(bigint *bi, size_t bit, int value);
size_t bi_used(const bigint *bi);
bi_leaf bi_to_intu(const bigint *bi);
bi_signed_leaf bi_to_int(const bigint *bi);
bi_intmax bi_to_intl(const bigint *bi);
bi_uintmax bi_to_intlu(const bigint *bi);

#ifdef BIGINT_ENABLE_LIBMATH
#include <float.h>
#include <math.h>

typedef enum
{
    bi_libmath_none,
    bi_libmath_inf,
    bi_libmath_nan
} bi_libmath_err;

bigint *bi_assign_float(bigint *bi, float f, bi_libmath_err *err);
bigint *bi_assign_double(bigint *bi, double f, bi_libmath_err *err);
bigint *bi_assign_doublel(bigint *bi, long double f, bi_libmath_err *err);
float bi_to_float(const bigint *bi);
double bi_to_double(const bigint *bi);
long double bi_to_doublel(const bigint *bi);
#endif // BIGINT_ENABLE_LIBMATH

typedef enum
{
    bi_rand_stdrand
} bi_rand_source;

bigint *bi_randgen(size_t bits, bi_rand_source source);
bigint *bi_randgen_range(const bigint *min, const bigint *max, bi_rand_source source);

bigint *bi_clear(bigint *bi);
int bi_is_negative(const bigint *bi);
void bi_set_negative(bigint *bi, int value);
int bi_is_growable(const bigint *bi);
void bi_set_growable(bigint *bi, int value);
int bi_is_freeable(const bigint *bi);
void bi_set_freeable(bigint *bi, int value);
int bi_is_destroyable(const bigint *bi);
void bi_set_destroyable(bigint *bi, int value);
int bi_cmp_mag(const bigint *bi, const bigint *bi2);
int bi_cmp(const bigint *bi, const bigint *bi2);
int bi_cmp_imm(const bigint *bi, bi_signed_leaf val);
int bi_cmp_immu(const bigint *bi, bi_leaf val);
int bi_cmp_zero(const bigint *bi);
int bi_is_zero(const bigint *bi);
int bi_is_one(const bigint *bi);
int bi_is_power_of_2(const bigint *bi);
int bi_log2(const bigint *bi);
#ifdef BIGINT_ENABLE_LIBMATH
long long bi_log2l(const bigint *bi);
int bi_log10_approx(const bigint *bi);
long long bi_log10l_approx(const bigint *bi);
int bi_log10(const bigint *bi);
long long bi_log10l(const bigint *bi);
int bi_logn_approx(const bigint *bi, uintmax_t n);
long long bi_lognl_approx(const bigint *bi, uintmax_t n);
int bi_logn(const bigint *bi, uintmax_t n);
long long bi_lognl(const bigint *bi, uintmax_t n);
#endif
int bi_trailing_zeroes(const bigint *bi);
bigint *bi_negate(const bigint *bi);
bigint *bi_negate_assign(bigint *bi);
bigint *bi_abs(const bigint *bi);
bigint *bi_abs_assign(bigint *bi);
bigint *bi_shl_leaf(const bigint *bi, size_t n);
bigint *bi_shl_leaf_assign(bigint *bi, size_t n);
bigint *bi_shr_leaf(const bigint *bi, size_t n);
bigint *bi_shr_leaf_assign(bigint *bi, size_t n);
bigint *bi_shl(const bigint *bi, size_t n);
bigint *bi_shl_assign(bigint *bi, size_t n);
bigint *bi_shl1(const bigint *bi);
bigint *bi_shl1_assign(bigint *bi);
bigint *bi_shr(const bigint *bi, size_t n);
bigint *bi_shr_assign(bigint *bi, size_t n);
bigint *bi_and(const bigint *bi, const bigint *bi2);
bigint *bi_and_assign(bigint *bi, const bigint *bi2);
bigint *bi_or(const bigint *bi, const bigint *bi2);
bigint *bi_or_assign(bigint *bi, const bigint *bi2);
bigint *bi_xor(const bigint *bi, const bigint *bi2);
bigint *bi_xor_assign(bigint *bi, const bigint *bi2);
bigint *bi_add(const bigint *bi, const bigint *bi2);
bigint *bi_add_assign(bigint *bi, const bigint *bi2);
bigint *bi_sub(const bigint *bi, const bigint *bi2);
bigint *bi_sub_assign(bigint *bi, const bigint *bi2);
bigint *bi_add_immediate(const bigint *bi, bi_signed_leaf val);
bigint *bi_add_immediate_assign(bigint *bi, bi_signed_leaf val);
bigint *bi_sub_immediate(const bigint *bi, bi_signed_leaf val);
bigint *bi_sub_immediate_assign(bigint *bi, bi_signed_leaf val);
bigint *bi_mul(const bigint *bi, const bigint *bi2);
bigint *bi_mul_assign(bigint *bi, const bigint *bi2);
bigint *bi_mul_immediateu(const bigint *bi, bi_leaf val);
bigint *bi_mul_immediateu_assign(bigint *bi, bi_leaf val);
bigint *bi_mul_immediate(const bigint *bi, bi_signed_leaf val);
bigint *bi_mul_immediate_assign(bigint *bi, bi_signed_leaf val);
bigint *bi_square(const bigint *bi);
bigint *bi_square_assign(bigint *bi);
bigint *bi_sqrt(const bigint *bi);
bigint *bi_sqrt_assign(bigint *bi);
bigint *bi_fact(bi_uintmax n);
bigint *bi_fibonacci(bi_uintmax n);
bigint *bi_fibonacci_mod(const bigint *n, const bigint *mod);
int bi_is_fibonacci(const bigint *bi);
bigint *bi_lucas(bi_uintmax n);
bigint *bi_lucas_mod(const bigint *n, const bigint *mod);
int bi_is_prime_low_divisibility_test(const bigint *bi);
int bi_is_prime_naive(const bigint *bi);
int bi_is_prime_fermat(const bigint *bi, size_t base);
int bi_is_prime_fibonacci(const bigint *bi);
int bi_is_prime_lucas(const bigint *bi);
int bi_is_prime_selfridge(const bigint *bi);
int bi_is_prime_miller_rabin(const bigint *bi, size_t k, bi_rand_source source);
bigint *bi_uexp(const bigint *bi, bi_uintmax n);
bigint *bi_uexp_assign(bigint *bi, bi_uintmax n);
bigint *bi_exp(const bigint *bi, bi_intmax n);
bigint *bi_exp_assign(bigint *bi, bi_intmax n);
bigint *bi_uexp_mod(const bigint *bi, bi_uintmax n, const bigint *mod);
bigint *bi_uexp_mod_assign(bigint *bi, bi_uintmax n, const bigint *mod);
bigint *bi_exp_mod(const bigint *bi, bi_intmax n, const bigint *mod);
bigint *bi_exp_mod_assign(bigint *bi, bi_intmax n, const bigint *mod);
bigint *bi_large_exp_mod(const bigint *bi, const bigint *n, const bigint *mod);
bigint *bi_large_exp_mod_assign(bigint *bi, const bigint *n, const bigint *mod);
bigint *bi_divmod(const bigint *bi, const bigint *bi2, bigint **q, bigint **r);
bigint *bi_divmod_assign(bigint *bi, const bigint *bi2, bigint **r);
bigint *bi_div(const bigint *bi, const bigint *bi2);
bigint *bi_div_assign(bigint *bi, const bigint *bi2);
bigint *bi_mod(const bigint *bi, const bigint *bi2);
bigint *bi_mod_assign(bigint *bi, const bigint *bi2);
bigint *bi_divmod_3(const bigint *bi, bigint **q, bi_signed_leaf *r);
bigint *bi_divmod_3_assign(bigint *bi, bi_signed_leaf *r);
bigint *bi_div_3(const bigint *bi);
bigint *bi_div_3_assign(bigint *bi);
bi_signed_leaf bi_mod_3(const bigint *bi);
bigint *bi_mod_3_assign(bigint *bi);
bigint *bi_divmod_10(const bigint *bi, bigint **q, bi_signed_leaf *r);
bigint *bi_divmod_10_assign(bigint *bi, bi_signed_leaf *r);
bigint *bi_div_10(const bigint *bi);
bigint *bi_div_10_assign(bigint *bi);
bi_signed_leaf bi_mod_10(const bigint *bi);
bigint *bi_mod_10_assign(bigint *bi);
bigint *bi_divmod_immediate(const bigint *bi, bi_signed_leaf denom, bigint **q, bi_signed_leaf *r);
bigint *bi_divmod_immediate_assign(bigint *bi, bi_signed_leaf denom, bi_signed_leaf *r);
bigint *bi_div_immediate(const bigint *bi, bi_signed_leaf denom);
bigint *bi_div_immediate_assign(bigint *bi, bi_signed_leaf denom);
bi_signed_leaf bi_mod_immediate(const bigint *bi, bi_signed_leaf denom, int *success);
bigint *bi_mod_immediate_assign(bigint *bi, bi_signed_leaf denom);
bigint *bi_gcd(const bigint *bi, const bigint *bi2);
void bi_swap(bigint *bi_a, bigint *bi_b);
int bi_sscan(const char *str, bigint *bi, size_t base);
int bi_sscan_n(const char *str, size_t len, bigint *bi, size_t base);
int bi_fscan(FILE *f, bigint *bi, size_t base);
int bi_scan(bigint *bi, size_t base);
bigint_string bi_sprint(const bigint *bi, size_t base);
void bi_fprint(FILE *f, const bigint *bi, size_t base);
void bi_print(const bigint *bi, size_t base);
void bi_destroy(bigint *bi);

#ifdef __cplusplus
} // extern "C"

#include <string>
#include <iostream>

// This class is safely destroyable if out of memory, but should not be operated on any further
class Bigint
{
    friend class Bigfrac;

    bigint d;

    Bigint(bigint *d) : d(*d)
    {
        this->d.flags &= ~BIGINT_FLAG_DESTROYABLE; // Make local copy indestructible so we don't free unallocated memory
        d->flags &= ~BIGINT_FLAG_FREEABLE; // Make data not freeable in sent copy, so we can keep a reference in the local copy without it being destroyed
        bi_destroy(d);
    }

public:
    struct error {};
    struct out_of_memory: public error {};
    struct division_by_zero: public error {};
    struct prime_gen_failed: public error {};

    const bigint *native_handle() const {return &d;}
    bigint *native_handle() {return &d;}

    Bigint() {if (bi_init(&d) == NULL) throw out_of_memory();}
    Bigint(int n) {if (bi_init_valuel(&d, n) == NULL) throw out_of_memory();}
    Bigint(bi_intmax n) {if (bi_init_valuel(&d, n) == NULL) throw out_of_memory();}
    Bigint(const char *s) {if (bi_init(&d) == NULL || bi_sscan(s, &d, 10) < 0) throw out_of_memory();}
    Bigint(const std::string &s) {if (bi_init(&d) == NULL || bi_sscan_n(s.c_str(), s.size(), &d, 10) < 0) throw out_of_memory();}
    Bigint(const Bigint &other) {if (bi_init_copy(&d, &other.d) == NULL) throw out_of_memory();}
#if __cplusplus >= 201103L
    Bigint(Bigint &&other) : d(std::move(other.d)) {other.d.data = NULL;}
#endif
    ~Bigint() {bi_destroy(&d);}

    Bigint &operator=(const Bigint &other)
    {
        if (bi_copy_to(&d, &other.d) == NULL) throw out_of_memory();
        return *this;
    }

    template<typename T>
    Bigint &operator=(T n)
    {
        if (bi_assignl(&d, n) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &operator=(bi_intmax n)
    {
        if (bi_assignl(&d, n) == NULL) throw out_of_memory();
        return *this;
    }

    void clear() {bi_clear(&d);}

    bool isNegative() const {return bi_is_negative(&d);}
    bool isZero() const {return bi_is_zero(&d);}
    bool isOne() const {return bi_is_one(&d);}
    bool isPowerOf2() const {return bi_is_power_of_2(&d);}
    // returns floor(log2(bigint)), -1 if input is zero
    int log2() const {return bi_log2(&d);}
    // returns floor(log2(bigint)), -1 if input is zero
    long long log2l() const {return bi_log2l(&d);}
#ifdef BIGINT_ENABLE_LIBMATH
    struct floating_point_error: public error {};
    struct floating_point_infinite: public floating_point_error {};
    struct floating_point_nan: public floating_point_error {};

    Bigint &operator=(float f)
    {
        bi_libmath_err err;

        if (bi_assign_float(&d, f, &err) == NULL) throw out_of_memory();

        if (err == bi_libmath_inf) throw floating_point_infinite();
        else if (err == bi_libmath_nan) throw floating_point_nan();

        return *this;
    }

    Bigint &operator=(double f)
    {
        bi_libmath_err err;

        if (bi_assign_double(&d, f, &err) == NULL) throw out_of_memory();

        if (err == bi_libmath_inf) throw floating_point_infinite();
        else if (err == bi_libmath_nan) throw floating_point_nan();

        return *this;
    }

    Bigint &operator=(long double f)
    {
        bi_libmath_err err;

        if (bi_assign_doublel(&d, f, &err) == NULL) throw out_of_memory();

        if (err == bi_libmath_inf) throw floating_point_infinite();
        else if (err == bi_libmath_nan) throw floating_point_nan();

        return *this;
    }

    // returns this value converted to float, returns +-INFINITY if out-of-range
    float toFloat() const
    {
        float f = bi_to_float(&d);
        if (isnan(f)) throw out_of_memory();
        return f;
    }
    // returns this value converted to double, returns +-INFINITY if out-of-range
    double toDouble() const
    {
        double f = bi_to_double(&d);
        if (isnan(f)) throw out_of_memory();
        return f;
    }
    // returns this value converted to long double, returns +-INFINITY if out-of-range
    long double toLDouble() const
    {
        long double f = bi_to_doublel(&d);
        if (isnan(f)) throw out_of_memory();
        return f;
    }
    // returns this value converted to an intmax_t, undefined if out-of-range
    intmax_t toInt() const {return bi_to_intl(&d);}
    // returns this value (disregarding sign) converted to a uintmax_t, undefined if out-of-range
    uintmax_t toUInt() const {return bi_to_intlu(&d);}

    // returns floor(log10(bigint)), -1 if input is zero
    int log10() const
    {
        int r = bi_log10(&d);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns floor(log10(bigint)), -1 if input is zero
    long long log10l() const
    {
        long long r = bi_log10l(&d);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns approximation of floor(log10(bigint)), -1 if input is zero
    // equal to floor(log10(bigint)) or floor(log10(bigint))-1
    int log10_approx() const
    {
        int r = bi_log10_approx(&d);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns approximation of floor(log10(bigint)), -1 if input is zero
    // equal to floor(log10(bigint)) or floor(log10(bigint))-1
    long long log10l_approx() const
    {
        long long r = bi_log10l_approx(&d);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns floor(logn(bigint)), -1 if input is zero
    int logn(uintmax_t n) const
    {
        int r = bi_logn(&d, n);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns floor(logn(bigint)), -1 if input is zero
    long long lognl(uintmax_t n) const
    {
        long long r = bi_lognl(&d, n);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns approximation of floor(logn(bigint)), -1 if input is zero
    // equal to floor(logn(bigint)) or floor(logn(bigint))-1
    int logn_approx(uintmax_t n) const
    {
        int r = bi_logn_approx(&d, n);
        if (r < -1) throw out_of_memory();
        return r;
    }
    // returns approximation of floor(logn(bigint)), -1 if input is zero
    // equal to floor(logn(bigint)) or floor(logn(bigint))-1
    long long lognl_approx(uintmax_t n) const
    {
        long long r = bi_lognl_approx(&d, n);
        if (r < -1) throw out_of_memory();
        return r;
    }
#endif
    // returns -1 if input is zero
    int trailingZeroes() const {return bi_trailing_zeroes(&d);}
    int compare(const Bigint &other) const {return bi_cmp(&d, &other.d);}
    int compareAbs(const Bigint &other) const {return bi_cmp_mag(&d, &other.d);}
    int signum() const {return bi_cmp_zero(&d);}

    Bigint operator+() const {return *this;}

    Bigint &negate() {bi_negate_assign(&d); return *this;}
    Bigint negated() const {return Bigint(*this).negate();}
    Bigint operator-() const {return Bigint(*this).negate();}

    Bigint &removeSign() {bi_abs_assign(&d); return *this;}
    Bigint abs() const {return Bigint(*this).removeSign();}

    // post- operators
    Bigint operator++(int)
    {
        Bigint b(*this);
        if (bi_add_immediate_assign(&d, 1) == NULL) throw out_of_memory();
        return b;
    }
    Bigint operator--(int)
    {
        Bigint b(*this);
        if (bi_sub_immediate_assign(&d, 1) == NULL) throw out_of_memory();
        return b;
    }

    // pre- operators
    Bigint &operator++()
    {
        if (bi_add_immediate_assign(&d, 1) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint &operator--()
    {
        if (bi_sub_immediate_assign(&d, 1) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &shiftLeft(size_t bits)
    {
        if (bi_shl_assign(&d, bits) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint shiftedLeft(size_t bits) const {return Bigint(*this).shiftLeft(bits);}
    Bigint &operator<<=(size_t bits) {return shiftLeft(bits);}
    Bigint operator<<(size_t bits) const {return Bigint(*this).shiftLeft(bits);}

    Bigint &shiftRight(size_t bits)
    {
        if (bi_shr_assign(&d, bits) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint shiftedRight(size_t bits) const {return Bigint(*this).shiftRight(bits);}
    Bigint &operator>>=(size_t bits) {return shiftRight(bits);}
    Bigint operator>>(size_t bits) const {return Bigint(*this).shiftRight(bits);}

    Bigint &add(const Bigint &other)
    {
        if (bi_add_assign(&d, &other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint added(const Bigint &other) const {return Bigint(*this).add(other);}
    Bigint &operator+=(const Bigint &other) {return add(other);}
    friend Bigint operator+(const Bigint &lhs, const Bigint &rhs) {return lhs.added(rhs);}

    Bigint &subtract(const Bigint &other)
    {
        if (bi_sub_assign(&d, &other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint subtracted(const Bigint &other) const {return Bigint(*this).subtract(other);}
    Bigint &operator-=(const Bigint &other) {return subtract(other);}
    friend Bigint operator-(const Bigint &lhs, const Bigint &rhs) {return lhs.subtracted(rhs);}

    Bigint &multiplyBy(const Bigint &other)
    {
        if (bi_mul_assign(&d, &other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint multipliedBy(const Bigint &other) const {return Bigint(*this).multiplyBy(other);}
    Bigint &operator*=(const Bigint &other) {return multiplyBy(other);}
    friend Bigint operator*(const Bigint &lhs, const Bigint &rhs) {return lhs.multipliedBy(rhs);}

    Bigint &square()
    {
        if (bi_square_assign(&d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint squared() const {return Bigint(*this).square();}

    Bigint &squareRoot()
    {
        if (bi_sqrt_assign(&d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint squareRooted() const {return Bigint(*this).squareRoot();}

    static Bigint rand(size_t bits, bi_rand_source source)
    {
        bigint *d = bi_randgen(bits, source);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    static Bigint rand(const Bigint &min, const Bigint &max, bi_rand_source source)
    {
        bigint *d = bi_randgen_range(&min.d, &max.d, source);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    static Bigint generateProbablePrime(size_t bits, size_t k, bi_rand_source source)
    {
        size_t tries = 100 * Bigint((bi_intmax) bits).log2();

        while (tries-- > 0)
        {
            bigint *d = bi_randgen(bits, source);
            if (d == NULL) throw out_of_memory();

            d->data[0] |= 1;

            if (bi_cmp_immu(d, 3) > 0)
            {
                int result = bi_is_prime_miller_rabin(d, k, source);
                if (result < 0)
                    throw out_of_memory();
                else if (result == 1)
                    return Bigint(d);
            }

            bi_destroy(d);
        }

        throw prime_gen_failed();
    }

    static Bigint factorial(bi_uintmax n)
    {
        bigint *d = bi_fact(n);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    static Bigint fibonacci(bi_uintmax n)
    {
        bigint *d = bi_fibonacci(n);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    static Bigint fibonacciMod(const Bigint &fibonacci, const Bigint &mod)
    {
        bigint *d = bi_fibonacci_mod(&fibonacci.d, &mod.d);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    static bool isFibonacci(const Bigint &n)
    {
        int r = bi_is_fibonacci(&n.d);
        if (r < 0) throw out_of_memory();
        return r;
    }

    static Bigint lucas(bi_uintmax n)
    {
        bigint *d = bi_lucas(n);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    static Bigint lucasMod(const Bigint &lucas, const Bigint &mod)
    {
        bigint *d = bi_lucas_mod(&lucas.d, &mod.d);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    bool isPrimeNaive() const
    {
        int r = bi_is_prime_naive(&d);
        if (r < 0) throw out_of_memory();
        return r;
    }

    bool isPrimeFermat(size_t base = 2) const
    {
        int r = bi_is_prime_fermat(&d, base);
        if (r < 0) throw out_of_memory();
        return r;
    }

    bool isPrimeFibonacci() const
    {
        int r = bi_is_prime_fibonacci(&d);
        if (r < 0) throw out_of_memory();
        return r;
    }

    bool isPrimeLucas() const
    {
        int r = bi_is_prime_lucas(&d);
        if (r < 0) throw out_of_memory();
        return r;
    }

    bool isPrimeSelfridge() const
    {
        int r = bi_is_prime_selfridge(&d);
        if (r < 0) throw out_of_memory();
        return r;
    }

    Bigint &powerU(bi_uintmax n)
    {
        if (bi_uexp_assign(&d, n) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint &power(bi_intmax n)
    {
        if (bi_exp_assign(&d, n) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint poweredU(bi_uintmax n) const {return Bigint(*this).power(n);}
    Bigint powered(bi_intmax n) const {return Bigint(*this).power(n);}
    Bigint &powerModU(bi_uintmax n, const Bigint &mod)
    {
        if (bi_uexp_mod_assign(&d, n, &mod.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint &powerMod(bi_intmax n, const Bigint &mod)
    {
        if (bi_exp_mod_assign(&d, n, &mod.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint &powerMod(const Bigint &n, const Bigint &mod)
    {
        if (bi_large_exp_mod_assign(&d, &n.d, &mod.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint poweredModU(bi_uintmax n, const Bigint &mod) const {return Bigint(*this).powerModU(n, mod);}
    Bigint poweredMod(bi_intmax n, const Bigint &mod) const {return Bigint(*this).powerMod(n, mod);}
    Bigint poweredMod(const Bigint &n, const Bigint &mod) const {return Bigint(*this).powerMod(n, mod);}

    /* TODO: Newton/Raphson division is currently very slow */
    Bigint newtonRaphson(const Bigint &other)
    {
        size_t k = log2()+1 + other.log2()+1;

        Bigint x(other.log2()+1), k2(2);

        k2 <<= k;

#if 0
        Bigint last, lastlast;
        while (true)
        {
            x *= k2 - x * other;
            x >>= k;
            if (x == last || x == lastlast)
                break;
            lastlast.swap(last);
            last = x;
        }
#else
        for (size_t i = 0; i < size_t(other.log2()+1); ++i)
            x = (x * (k2 - x * other)) >> k;
#endif

        x *= *this;
        x >>= k;

        if (*this - (x * other) >= other)
            ++x;

        return x;
    }

    Bigint &divideBy(const Bigint &other)
    {
        if (bi_is_zero(&other.d)) throw division_by_zero();
        if (bi_div_assign(&d, &other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint divided(const Bigint &other) const {return Bigint(*this).divideBy(other);}
    Bigint &operator/=(const Bigint &other) {return divideBy(other);}
    friend Bigint operator/(const Bigint &lhs, const Bigint &rhs) {return lhs.divided(rhs);}

    Bigint &moduloBy(const Bigint &other)
    {
        if (bi_is_zero(&other.d)) throw division_by_zero();
        if (bi_mod_assign(&d, &other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint modulo(const Bigint &other) const {return Bigint(*this).moduloBy(other);}
    Bigint &operator%=(const Bigint &other) {return moduloBy(other);}
    friend Bigint operator%(const Bigint &lhs, const Bigint &rhs) {return lhs.modulo(rhs);}

    Bigint &gcdAssign(const Bigint &other)
    {
        bigint *n = bi_gcd(&d, &other.d);
        if (n == NULL) throw out_of_memory();
        bi_swap(&d, n);
        bi_destroy(n);
        return *this;
    }
    Bigint gcd(const Bigint &other) const {return Bigint(*this).gcdAssign(other);}

    static Bigint fromString(const char *s, size_t base = 10, int *charsRead = NULL)
    {
        bigint *b = bi_new();
        int chars;
        if (b == NULL ||
            (chars = bi_sscan(s, b, base)) < 0) throw out_of_memory();

        if (charsRead != NULL)
            *charsRead = chars;

        return Bigint(b);
    }

    static Bigint fromString(const std::string &s, size_t base = 10, int *charsRead = NULL)
    {
        bigint *b = bi_new();
        int chars;
        if (b == NULL ||
            (chars = bi_sscan_n(s.c_str(), s.size(), b, base)) < 0) throw out_of_memory();

        if (charsRead != NULL)
            *charsRead = chars;

        return Bigint(b);
    }

    std::string toString(size_t base = 10) const
    {
        bigint_string bstr = bi_sprint(&d, base);
        if (bstr.string == NULL) throw out_of_memory();

        std::string str(bstr.string, bstr.len);
        bis_destroy(bstr);
        return str;
    }

    void swap(Bigint &other) {bi_swap(&d, &other.d);}

    friend bool operator<(const Bigint &lhs, const Bigint &rhs) {return lhs.compare(rhs) < 0;}
    friend bool operator>(const Bigint &lhs, const Bigint &rhs) {return lhs.compare(rhs) > 0;}
    friend bool operator<=(const Bigint &lhs, const Bigint &rhs) {return lhs.compare(rhs) <= 0;}
    friend bool operator>=(const Bigint &lhs, const Bigint &rhs) {return lhs.compare(rhs) >= 0;}
    friend bool operator==(const Bigint &lhs, const Bigint &rhs) {return lhs.compare(rhs) == 0;}
    friend bool operator!=(const Bigint &lhs, const Bigint &rhs) {return lhs.compare(rhs) != 0;}

    friend std::istream &operator>>(std::istream &in, Bigint &bi)
    {
        std::ios_base::fmtflags flags = in.flags();
        std::string s;
        in >> s;

        switch (flags & std::ios_base::basefield)
        {
            case std::ios_base::hex:
                if (s.size() >= 2 && s[0] == '0' && tolower(s[1]) == 'x')
                    s.erase(0, 2);
                bi = fromString(s, 16);
                return in;
            case std::ios_base::oct:
                bi = fromString(s, 8);
                return in;
            case std::ios_base::dec:
            default:
                bi = fromString(s, 10);
                return in;
        }
    }

    friend std::ostream &operator<<(std::ostream &out, const Bigint &bi)
    {
        std::ios_base::fmtflags flags = out.flags();

        switch (flags & std::ios_base::basefield)
        {
            case std::ios_base::hex:
                return out << "0x" << bi.toString(16);
            case std::ios_base::oct:
            {
                std::string str = bi.toString(8);
                if (str == "0")
                    return out << str;
                else
                    return out << '0' << str;
            }
            case std::ios_base::dec:
            default:
                return out << bi.toString(10);
        }
    }
};

inline void swap(Bigint &me, Bigint &other)
{
    me.swap(other);
}

#endif

#endif // BIGINT_H
