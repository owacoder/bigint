#ifndef BIGINT_H
#define BIGINT_H

#include <memory.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>

#if 0// defined(__GNUC__)
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
#else
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
#endif
// Minimum number of leaves allocated for a bigint. Must be greater than 0.
#define BIGINT_MINLEAFS (8)

#define BIGINT_LEAFMAX ((bi_leaf) ~0ULL)
#define BIGINT_SIGNED_LEAFMAX ((bi_signed_leaf) (BIGINT_LEAFMAX >> 1))
#define BIGINT_LEAFBYTES (BIGINT_LEAFSIZE*CHAR_BIT/8)
#define BIGINT_LEAFBITS (BIGINT_LEAFSIZE*CHAR_BIT)
#define BIGINT_LEAFS_PER_BI_INTMAX (sizeof(bi_intmax)/BIGINT_LEAFSIZE)
#define BIGINT_KARATSUBA_MIN_LEAFS (640/BIGINT_LEAFBITS)
#define BIGINT_TOOM_COOK_MIN_LEAFS (10000/BIGINT_LEAFBITS)

#define BIGINT_ENABLE_PTHREADS
#define BIGINT_ENABLE_LIBMATH

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

typedef struct
{
    bi_leaf *data;
    size_t size, sign; //sign: 0 = positive, 1 = negative
} bigint;

bigint *bi_new();
bigint *bi_copy(const bigint *bi);
bigint *bi_copy_mag_to(bigint *dst, const bigint *src);
bigint *bi_copy_to(bigint *dst, const bigint *src);
bigint *bi_assignu(bigint *bi, bi_leaf value);
bigint *bi_assign(bigint *bi, bi_signed_leaf value);
bigint *bi_assignl(bigint *bi, bi_intmax value);
bigint *bi_assignlu(bigint *bi, bi_uintmax value);
size_t bi_used(const bigint *bi);
bi_leaf bi_to_intu(const bigint *bi);
bi_signed_leaf bi_to_int(const bigint *bi);
bi_intmax bi_to_intl(const bigint *bi);
bi_uintmax bi_to_intlu(const bigint *bi);
bigint *bi_clear(bigint *bi);
int bi_is_negative(const bigint *bi);
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
int bi_is_power_of_10(const bigint *bi);
int bi_log10_approx(const bigint *bi);
int bi_log10(const bigint *bi);
int bi_logn_approx(const bigint *bi, uintmax_t n);
int bi_logn(const bigint *bi, uintmax_t n);
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
bigint *bi_uexp(const bigint *bi, bi_uintmax n);
bigint *bi_uexp_assign(bigint *bi, bi_uintmax n);
bigint *bi_exp(const bigint *bi, bi_intmax n);
bigint *bi_exp_assign(bigint *bi, bi_intmax n);
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
bi_signed_leaf bi_mod_immediate(const bigint *bi, bi_signed_leaf denom);
bigint *bi_mod_immediate_assign(bigint *bi, bi_signed_leaf denom);
bigint *bi_gcd(const bigint *bi, const bigint *bi2);
void bi_swap(bigint *bi_a, bigint *bi_b);
int bi_sscan(const char *str, bigint *bi, size_t base);
int bi_fscan(FILE *f, bigint *bi, size_t base);
int bi_scan(bigint *bi, size_t base);
bigint_string bi_sprint(const bigint *bi, size_t base);
void bi_fprint(FILE *f, const bigint *bi, size_t base);
void bi_print(const bigint *bi, size_t base);
void bi_destroy(bigint *bi);

#ifdef __cplusplus
} // extern "C"

#include <string>

class Bigint
{
    bigint *d;

    Bigint(bigint *d) : d(d) {}

public:
    struct error {};
    struct out_of_memory: public error {};
    struct division_by_zero: public error {};

    const bigint *native_handle() const {return d;}
    bigint *&native_handle() {return d;}

    Bigint() : d(bi_new()) {if (d == NULL) throw out_of_memory();}
    Bigint(bi_leaf n) : d(bi_new()) {if (bi_assignu(d, n) == NULL) throw out_of_memory();}
    Bigint(bi_signed_leaf n) : d(bi_new()) {if (bi_assign(d, n) == NULL) throw out_of_memory();}
    Bigint(bi_intmax n) : d(bi_new()) {if (bi_assignl(d, n) == NULL) throw out_of_memory();}
    Bigint(bi_uintmax n) : d(bi_new()) {if (bi_assignlu(d, n) == NULL) throw out_of_memory();}
    Bigint(const Bigint &other) : d(bi_copy(other.d)) {if (d == NULL) throw out_of_memory();}
#if __cplusplus >= 201103L
    Bigint(Bigint &&other) : d(other.d) {other.d = NULL;}
#endif
    ~Bigint() {bi_destroy(d);}

    Bigint &operator=(const Bigint &other)
    {
        if (&other == this) return *this;
        if (bi_copy_to(d, other.d) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &operator=(bi_leaf n)
    {
        if (bi_assignu(d, n) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &operator=(bi_signed_leaf n)
    {
        if (bi_assign(d, n) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &operator=(bi_intmax n)
    {
        if (bi_assignl(d, n) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &operator=(bi_uintmax n)
    {
        if (bi_assignlu(d, n) == NULL) throw out_of_memory();
        return *this;
    }

    void clear() {bi_clear(d);}

    bool isNegative() const {return bi_is_negative(d);}
    bool isZero() const {return bi_is_zero(d);}
    bool isOne() const {return bi_is_one(d);}
    bool isPowerOf2() const {return bi_is_power_of_2(d);}
    // returns floor(log2(bigint)), -1 if input is zero
    int log2() const {return bi_log2(d);}
#ifdef BIGINT_ENABLE_LIBMATH
    // returns floor(log10(bigint)), -1 if input is zero
    int log10() const
    {
        int r = bi_log10(d);
        return r < 0? -1: r;
    }
    int logn(uintmax_t n) const
    {
        int r = bi_logn(d, n);
        return r < 0? -1: r;
    }
#endif
    // returns -1 if input is zero
    int trailingZeroes() const {return bi_trailing_zeroes(d);}
    int compare(const Bigint &other) const {return bi_cmp(d, other.d);}
    int compareAbs(const Bigint &other) const {return bi_cmp_mag(d, other.d);}
    int signum() const {return bi_cmp_zero(d);}

    Bigint operator+() const {return *this;}

    Bigint &negate() {bi_negate(d); return *this;}
    Bigint negated() const {return Bigint(*this).negate();}
    Bigint operator-() const {return Bigint(*this).negate();}

    Bigint &removeSign() {bi_abs_assign(d); return *this;}
    Bigint abs() const {return Bigint(*this).removeSign();}

    // post- operators
    Bigint operator++(int)
    {
        Bigint b(*this);
        if (bi_add_immediate_assign(d, 1) == NULL) throw out_of_memory();
        return b;
    }
    Bigint operator--(int)
    {
        Bigint b(*this);
        if (bi_sub_immediate_assign(d, 1) == NULL) throw out_of_memory();
        return b;
    }

    // pre- operators
    Bigint &operator++()
    {
        if (bi_add_immediate_assign(d, 1) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint &operator--()
    {
        if (bi_sub_immediate_assign(d, 1) == NULL) throw out_of_memory();
        return *this;
    }

    Bigint &shiftLeft(size_t bits)
    {
        if (bi_shl_assign(d, bits) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint shiftedLeft(size_t bits) const {return Bigint(*this).shiftLeft(bits);}
    Bigint &operator<<=(size_t bits) {return shiftLeft(bits);}
    Bigint operator<<(size_t bits) const {return Bigint(*this).shiftLeft(bits);}

    Bigint &shiftRight(size_t bits)
    {
        if (bi_shr_assign(d, bits) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint shiftedRight(size_t bits) const {return Bigint(*this).shiftRight(bits);}
    Bigint &operator>>=(size_t bits) {return shiftRight(bits);}
    Bigint operator>>(size_t bits) const {return Bigint(*this).shiftRight(bits);}

    Bigint &add(const Bigint &other)
    {
        if (bi_add_assign(d, other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint added(const Bigint &other) const {return Bigint(*this).add(other);}
    Bigint &operator+=(const Bigint &other) {return add(other);}
    Bigint operator+(const Bigint &other) const {return Bigint(*this).add(other);}

    Bigint &subtract(const Bigint &other)
    {
        if (bi_sub_assign(d, other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint subtracted(const Bigint &other) const {return Bigint(*this).subtract(other);}
    Bigint &operator-=(const Bigint &other) {return subtract(other);}
    Bigint operator-(const Bigint &other) const {return Bigint(*this).subtract(other);}

    Bigint &multiplyBy(const Bigint &other)
    {
        if (bi_mul_assign(d, other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint multiplied(const Bigint &other) const {return Bigint(*this).multiplyBy(other);}
    Bigint &operator*=(const Bigint &other) {return multiplyBy(other);}
    Bigint operator*(const Bigint &other) const {return Bigint(*this).multiplyBy(other);}

    Bigint &square()
    {
        if (bi_square_assign(d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint squared() const {return Bigint(*this).square();}

    Bigint &squareRoot()
    {
        if (bi_sqrt_assign(d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint squareRooted() const {return Bigint(*this).squareRoot();}

    static Bigint factorial(bi_uintmax n)
    {
        bigint *d = bi_fact(n);
        if (d == NULL) throw out_of_memory();
        return Bigint(d);
    }

    Bigint &power(bi_uintmax n)
    {
        if (bi_uexp_assign(d, n) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint &power(bi_intmax n)
    {
        if (bi_exp_assign(d, n) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint powered(bi_uintmax n) const {return Bigint(*this).power(n);}
    Bigint powered(bi_intmax n) const {return Bigint(*this).power(n);}

    Bigint &divideBy(const Bigint &other)
    {
        if (bi_is_zero(other.d)) throw division_by_zero();
        if (bi_div_assign(d, other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint divided(const Bigint &other) const {return Bigint(*this).divideBy(other);}
    Bigint &operator/=(const Bigint &other) {return divideBy(other);}
    Bigint operator/(const Bigint &other) const {return Bigint(*this).divideBy(other);}

    Bigint &moduloBy(const Bigint &other)
    {
        if (bi_is_zero(other.d)) throw division_by_zero();
        if (bi_mod_assign(d, other.d) == NULL) throw out_of_memory();
        return *this;
    }
    Bigint modulo(const Bigint &other) const {return Bigint(*this).moduloBy(other);}
    Bigint &operator%=(const Bigint &other) {return moduloBy(other);}
    Bigint operator%(const Bigint &other) const {return Bigint(*this).moduloBy(other);}

    Bigint &gcdAssign(const Bigint &other)
    {
        bigint *n = bi_gcd(d, other.d);
        if (n == NULL) throw out_of_memory();
        bi_swap(d, n);
        bi_destroy(n);
        return *this;
    }
    Bigint gcd(const Bigint &other) const {return Bigint(*this).gcdAssign(other);}

    std::string toString(size_t base = 10) const
    {
        bigint_string bstr = bi_sprint(d, base);
        if (bstr.string == NULL) throw out_of_memory();

        std::string str(bstr.string, bstr.len);
        bis_destroy(bstr);
        return str;
    }
};

#endif

#endif // BIGINT_H
