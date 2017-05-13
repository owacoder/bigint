#ifndef BIGFRAC_H
#define BIGFRAC_H

#include "bigint.h"

typedef struct
{
    bigint *n, *d;
} bigfrac;

#ifdef __cplusplus
extern "C" {
#endif

bigfrac *bf_new();
bigfrac *bf_new_value(bi_signed_leaf value);
bigfrac *bf_new_valueu(bi_leaf value);
bigfrac *bf_new_valuel(bi_intmax value);
bigfrac *bf_new_valuelu(bi_uintmax value);
bigfrac *bf_new_frac_value(bi_signed_leaf numer, bi_signed_leaf denom);
bigfrac *bf_new_frac_valueu(bi_leaf numer, bi_leaf denom);
bigfrac *bf_new_frac_valuel(bi_intmax numer, bi_intmax denom);
bigfrac *bf_new_frac_valuelu(bi_uintmax numer, bi_uintmax denom);
bigfrac *bf_new_frac(const bigint *numer, const bigint *denom);
bigfrac *bf_new_frac_init(bigint *numer, bigint *denom);
bigfrac *bf_copy(const bigfrac *bf);
bigfrac *bf_copy_to(bigfrac *dst, const bigfrac *src);
bigfrac *bf_copy_frac_to(bigfrac *dst, const bigint *numer, const bigint *denom);
bigfrac *bf_copy_frac_init_to(bigfrac *dst, bigint *numer, bigint *denom);
bigfrac *bf_copy_numer_to(bigfrac *bf, const bigint *numer);
bigfrac *bf_copy_denom_to(bigfrac *bf, const bigint *denom);
bigfrac *bf_clear(bigfrac *bf);
bigfrac *bf_assign(bigfrac *bf, bi_signed_leaf value);
bigfrac *bf_assignu(bigfrac *bf, bi_leaf value);
bigfrac *bf_assignl(bigfrac *bf, bi_intmax value);
bigfrac *bf_assignlu(bigfrac *bf, bi_uintmax value);
bigfrac *bf_frac_assign(bigfrac *bf, bi_signed_leaf numer, bi_signed_leaf denom);
bigfrac *bf_frac_assignu(bigfrac *bf, bi_leaf numer, bi_leaf denom);
bigfrac *bf_frac_assignl(bigfrac *bf, bi_intmax numer, bi_intmax denom);
bigfrac *bf_frac_assignlu(bigfrac *bf, bi_uintmax numer, bi_uintmax denom);
int bf_cmp_eq(const bigfrac *bf, const bigfrac *bf2);
int bf_cmp(const bigfrac *bf, const bigfrac *bf2);
int bf_cmp_zero(const bigfrac *bf);
int bf_is_negative(const bigfrac *bf);
int bf_is_one(const bigfrac *bf);
int bf_is_zero(const bigfrac *bf);
int bf_is_undefined(const bigfrac *bf);
int bf_is_power_of_2(const bigfrac *bf);
bigfrac *bf_negate(const bigfrac *bf);
bigfrac *bf_negate_assign(bigfrac *bf);
bigfrac *bf_reciprocal(const bigfrac *bf);
bigfrac *bf_reciprocal_assign(bigfrac *bf);
bigfrac *bf_abs(const bigfrac *bf);
bigfrac *bf_abs_assign(bigfrac *bf);
bigfrac *bf_add(const bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_add_assign(bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_add_immediate(const bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_add_immediate_assign(bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_sub(const bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_sub_assign(bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_sub_immediate(const bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_sub_immediate_assign(bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_mul(const bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_mul_assign(bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_mul_immediate(const bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_mul_immediate_assign(bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_mul_immediateu(const bigfrac *bf, bi_leaf val);
bigfrac *bf_mul_immediateu_assign(bigfrac *bf, bi_leaf val);
bigfrac *bf_square(const bigfrac *bf);
bigfrac *bf_square_assign(bigfrac *bf);
bigfrac *bf_sqrt(const bigfrac *bf);
bigfrac *bf_sqrt_assign(bigfrac *bf);
bigfrac *bf_div(const bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_div_assign(bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_div_immediate(const bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_div_immediate_assign(bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_div_immediateu(const bigfrac *bf, bi_leaf val);
bigfrac *bf_div_immediateu_assign(bigfrac *bf, bi_leaf val);
bigfrac *bf_reduce(bigfrac *bf);
int bf_sscan(const char *s, bigfrac *bf, size_t base);
int bf_sscan_n(const char *s, size_t len, bigfrac *bf, size_t base);
bigint_string bf_sprint(const bigfrac *bf, size_t base);
void bf_print(const bigfrac *bf, size_t base);
void bf_fprint(FILE *out, const bigfrac *bf, size_t base);
void bf_swap(bigfrac *bfa, bigfrac *bfb);
void bf_destroy(bigfrac *bf);

#ifdef __cplusplus
} // extern "C"

#include <string>
#include <iostream>

class Bigfrac
{
    bigfrac *d;

    Bigfrac(bigfrac *d) : d(d) {}
    Bigfrac(bigint *n, bigint *d) : d(bf_new_frac_init(n, d)) {if (d == NULL) throw out_of_memory();}

public:
    struct error {};
    struct out_of_memory: public error {};
    struct division_by_zero: public error {};

    const bigfrac *native_handle() const {return d;}
    bigfrac *&native_handle() {return d;}

    Bigfrac() : d(bf_new()) {if (d == NULL) throw out_of_memory();}
    Bigfrac(int n) : d(bf_new_valuel(n)) {if (d == NULL) throw out_of_memory();}
    Bigfrac(int numer, int denom) : d(bf_new_frac_valuel(numer, denom)) {if (d == NULL) throw out_of_memory();}
    Bigfrac(bi_intmax n) : d(bf_new_valuel(n)) {if (d == NULL) throw out_of_memory();}
    Bigfrac(bi_intmax numer, bi_intmax denom) : d(bf_new_frac_valuel(numer, denom)) {if (d == NULL) throw out_of_memory();}
    Bigfrac(const char *s) : d(bf_new()) {if (d == NULL || bf_sscan(s, d, 10) < 0) throw out_of_memory();}
    Bigfrac(const std::string &s) : d(bf_new()) {if (d == NULL || bf_sscan_n(s.c_str(), s.size(), d, 10) < 0) throw out_of_memory();}
    Bigfrac(const Bigint &other) : d(bf_new_frac_init(bi_copy(other.d), bi_new_valueu(1))) {if (d == NULL) throw out_of_memory();}
    Bigfrac(const Bigint &numer, const Bigint &denom) : d(bf_new_frac_init(bi_copy(numer.d), bi_copy(denom.d))) {if (d == NULL) throw out_of_memory();}
    Bigfrac(const Bigfrac &other) : d(bf_copy(other.d)) {if (d == NULL) throw out_of_memory();}
#if __cplusplus >= 201103L
    Bigfrac(Bigfrac &&other) : d(other.d) {other.d = NULL;}
#endif
    ~Bigfrac() {bf_destroy(d);}

    Bigfrac &operator=(const Bigfrac &other)
    {
        if ((d = bf_copy_to(d, other.d)) == NULL) throw out_of_memory();
        return *this;
    }

    Bigfrac &operator=(const Bigint &other)
    {
        if ((d->n = bi_copy_to(d->n, other.d)) == NULL ||
            (d->d = bi_assignu(d->d, 1)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
        return *this;
    }

    Bigfrac &operator=(bi_intmax n)
    {
        if ((d = bf_assignl(d, n)) == NULL) throw out_of_memory();
        return *this;
    }

    void setFraction(const Bigint &numer, const Bigint &denom)
    {
        setNumerator(numer);
        setDenominator(denom);
    }
    void setFraction(bi_intmax numer, bi_intmax denom)
    {
        setNumerator(numer);
        setDenominator(denom);
    }

    void setNumerator(const Bigint &numer)
    {
        if ((d->n = bi_copy_to(d->n, numer.d)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
    }
    void setNumerator(bi_intmax numer)
    {
        if ((d->n = bi_assignl(d->n, numer)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
    }
    void setDenominator(const Bigint &denom)
    {
        if ((d->d = bi_copy_to(d->d, denom.d)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
    }
    void setDenominator(bi_intmax denom)
    {
        if ((d->d = bi_assignl(d->d, denom)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
    }

    void clear() {bf_clear(d);}

    Bigint numerator() const
    {
        bigint *b = bi_copy(d->n);
        if (b == NULL) throw out_of_memory();
        return Bigint(b);
    }
    Bigint denominator() const
    {
        bigint *b = bi_copy(d->d);
        if (b == NULL) throw out_of_memory();
        return Bigint(b);
    }
    Bigint toBigint() const
    {
        bigint *b = bi_div(d->n, d->d);
        if (b == NULL) throw out_of_memory();
        return Bigint(b);
    }

    bool isNegative() const {return bf_is_negative(d);}
    bool isZero() const {return bf_is_zero(d);}
    bool isOne() const {return bf_is_one(d);}
    bool isUndefined() const {return bf_is_undefined(d);}
    bool isPowerOf2() const {return bf_is_power_of_2(d);}
    int compare(const Bigfrac &other) const
    {
        int n = bf_cmp(d, other.d);
        if (n == 2) throw out_of_memory();
        return n;
    }
    bool equalTo(const Bigfrac &other) const
    {
        int n = bf_cmp_eq(d, other.d);
        if (n == 2) throw out_of_memory();
        return n == 0;
    }
    int signum() const {return bf_cmp_zero(d);}

    Bigfrac operator+() const {return *this;}

    Bigfrac &negate() {bf_negate(d); return *this;}
    Bigfrac negated() const {return Bigfrac(*this).negate();}
    Bigfrac operator-() const {return Bigfrac(*this).negate();}

    Bigfrac &removeSign() {bf_abs_assign(d); return *this;}
    Bigfrac abs() const {return Bigfrac(*this).removeSign();}

    // post- operators
    Bigfrac operator++(int)
    {
        Bigfrac b(*this);
        if ((d = bf_add_immediate_assign(d, 1)) == NULL) throw out_of_memory();
        return b;
    }
    Bigfrac operator--(int)
    {
        Bigfrac b(*this);
        if ((d = bf_sub_immediate_assign(d, 1)) == NULL) throw out_of_memory();
        return b;
    }

    // pre- operators
    Bigfrac &operator++()
    {
        if ((d = bf_add_immediate_assign(d, 1)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac &operator--()
    {
        if ((d = bf_sub_immediate_assign(d, 1)) == NULL) throw out_of_memory();
        return *this;
    }

    Bigfrac &shiftLeft(size_t bits)
    {
        if ((d->n = bi_shl_assign(d->n, bits)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
        return *this;
    }
    Bigfrac shiftedLeft(size_t bits) const {return Bigfrac(*this).shiftLeft(bits);}
    Bigfrac &operator<<=(size_t bits) {return shiftLeft(bits);}
    Bigfrac operator<<(size_t bits) const {return Bigfrac(*this).shiftLeft(bits);}

    Bigfrac &shiftRight(size_t bits)
    {
        if ((d->n = bi_shr_assign(d->n, bits)) == NULL)
        {
            bf_destroy(d); d=NULL;
            throw out_of_memory();
        }
        return *this;
    }
    Bigfrac shiftedRight(size_t bits) const {return Bigfrac(*this).shiftRight(bits);}
    Bigfrac &operator>>=(size_t bits) {return shiftRight(bits);}
    Bigfrac operator>>(size_t bits) const {return Bigfrac(*this).shiftRight(bits);}

    Bigfrac &add(const Bigfrac &other)
    {
        if ((d = bf_add_assign(d, other.d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac added(const Bigfrac &other) const {return Bigfrac(*this).add(other);}
    Bigfrac &operator+=(const Bigfrac &other) {return add(other);}
    friend Bigfrac operator+(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.added(rhs);}

    Bigfrac &subtract(const Bigfrac &other)
    {
        if ((d = bf_sub_assign(d, other.d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac subtracted(const Bigfrac &other) const {return Bigfrac(*this).subtract(other);}
    Bigfrac &operator-=(const Bigfrac &other) {return subtract(other);}
    friend Bigfrac operator-(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.subtracted(rhs);}

    Bigfrac &multiplyBy(const Bigfrac &other)
    {
        if ((d = bf_mul_assign(d, other.d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac multipliedBy(const Bigfrac &other) const {return Bigfrac(*this).multiplyBy(other);}
    Bigfrac &operator*=(const Bigfrac &other) {return multiplyBy(other);}
    friend Bigfrac operator*(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.multipliedBy(rhs);}

    Bigfrac &square()
    {
        if ((d = bf_square_assign(d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac squared() const {return Bigfrac(*this).square();}

    Bigfrac &squareRoot()
    {
        if ((d = bf_sqrt_assign(d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac squareRooted() const {return Bigfrac(*this).squareRoot();}

    static Bigfrac factorial(bi_uintmax n)
    {
        bigint *d = bi_fact(n);
        if (d == NULL) throw out_of_memory();
        return Bigfrac(d, bi_new_valueu(1));
    }

    static Bigfrac fibonacci(bi_uintmax n)
    {
        bigint *d = bi_fibonacci(n);
        if (d == NULL) throw out_of_memory();
        return Bigfrac(d, bi_new_valueu(1));
    }

    static bool isFibonacci(const Bigfrac &n)
    {
        int r = bi_is_fibonacci(n.d->n);
        if (r < 0) throw out_of_memory();
        return r && bi_is_one(n.d->d);
    }

    Bigfrac &divideBy(const Bigfrac &other)
    {
        if (bf_is_zero(other.d)) throw division_by_zero();
        if ((d = bf_div_assign(d, other.d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac divided(const Bigfrac &other) const {return Bigfrac(*this).divideBy(other);}
    Bigfrac &operator/=(const Bigfrac &other) {return divideBy(other);}
    friend Bigfrac operator/(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.divided(rhs);}

    Bigfrac &reduce()
    {
        if ((d = bf_reduce(d)) == NULL) throw out_of_memory();
        return *this;
    }
    Bigfrac reduced() const {return Bigfrac(*this).reduce();}

    static Bigfrac fromString(const char *s, size_t base = 10)
    {
        bigfrac *b = bf_new();
        if (b == NULL ||
            bf_sscan(s, b, base) < 0) throw out_of_memory();
        return Bigfrac(b);
    }

    static Bigfrac fromString(const std::string &s, size_t base = 10)
    {
        bigfrac *b = bf_new();
        if (b == NULL ||
            bf_sscan_n(s.c_str(), s.size(), b, base) < 0) throw out_of_memory();
        return Bigfrac(b);
    }

    std::string toString(size_t base = 10) const
    {
        bigint_string bstr = bf_sprint(d, base);
        if (bstr.string == NULL) throw out_of_memory();

        std::string str(bstr.string, bstr.len);
        bis_destroy(bstr);
        return str;
    }

    friend bool operator<(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.compare(rhs) < 0;}
    friend bool operator>(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.compare(rhs) > 0;}
    friend bool operator<=(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.compare(rhs) <= 0;}
    friend bool operator>=(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.compare(rhs) >= 0;}
    friend bool operator==(const Bigfrac &lhs, const Bigfrac &rhs) {return lhs.equalTo(rhs);}
    friend bool operator!=(const Bigfrac &lhs, const Bigfrac &rhs) {return !lhs.equalTo(rhs);}

    friend std::istream &operator>>(std::istream &in, Bigfrac &bi)
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

    friend std::ostream &operator<<(std::ostream &out, const Bigfrac &bi)
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

#endif

#endif // BIGFRAC_H

