#include "bigint.h"
#include "general.h"
#include "io.h"

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#ifdef BIGINT_ENABLE_LIBMATH
#include <math.h>
#endif

#ifdef BIGINT_ENABLE_PTHREADS
#include <pthread.h>
#endif

#if BIGINT_MINLEAFS < 1
#error BIGINT_MINLEAFS must be at least 1
#endif

void mul128(uint64_t a, uint64_t b, uint64_t *h, uint64_t *l)
{
    /* final result is equal to ha*hb*2^64 + (la*hb+ha*lb)*2^32 + la*lb */
    uint64_t t, m;
    uint64_t la = a & 0xFFFFFFFF, lb = b & 0xFFFFFFFF, ha = a >> 32, hb = b >> 32;

    t = la * hb;
    *h = t >> 32;
    m = ha * lb;
    *h += m >> 32;
    *h += ha * hb;
    *l = a * b;
}

/* create a new bigint with given size initialized to 0 and return it */
/* returns NULL on out of memory condition */
static bigint *bi_new_sized(size_t size)
{
    bigint *ptr = malloc(sizeof(*ptr));
    if (ptr != NULL)
    {
        ptr->size = size;
        ptr->sign = 0;
        ptr->data = calloc(size, sizeof(*ptr->data));
        if (ptr->data == NULL)
        {
            free(ptr);
            return NULL;
        }
    }
    return ptr;
}

/* resize bi to given size, preserving the value, and return bi */
/* returns NULL and destroys bi on out of memory condition */
static bigint *bi_resize(bigint *bi, size_t size)
{
    void *temp;
    size_t old_size;
    if (bi->size >= size) return bi;

    old_size = bi->size;
    bi->size = size;
    temp = realloc(bi->data, size * sizeof(bi_leaf));
    if (temp == NULL)
    {
        free(bi->data);
        free(bi);
        return NULL;
    }
    bi->data = temp;
    memset(bi->data + old_size, 0, (size-old_size) * sizeof(bi_leaf));

    return bi;
}

/* returns 1 if the topmost leaf of bi is non-zero */
static size_t bi_full(const bigint *bi)
{
    return bi->size == 0 || bi->data[bi->size - 1] != 0;
}

/* returns bit position of highest set bit, 1-indexed
 * (equivalent to the number of bits required to store the number) */
static size_t bi_bitcount(const bigint *bi)
{
    size_t i = bi->size;

    while (i > 0)
        if (bi->data[--i] != 0)
        {
            bi_leaf l = bi->data[i];
            i *= BIGINT_LEAFBITS;
            do ++i; while (l >>= 1);
            return i;
        }

    return 0;
}

/* returns bit at specified position, 0-indexed */
/* returns -1 if specified bit is out of range */
bi_signed_leaf bi_bit(const bigint *bi, size_t bit)
{
    size_t word = bit / BIGINT_LEAFBITS;
    if (word >= bi->size)
        return -1;
    return (bi->data[word] >> (bit % BIGINT_LEAFBITS)) & 1;
}

/* returns number of leaves used in bi */
size_t bi_used(const bigint *bi)
{
    size_t i = bi->size;

    for (; i > 0; )
        if (bi->data[--i] != 0)
            return i+1;

    return 0;
}

/* create a new bigint initialized to 0 and return it */
/* returns NULL on out of memory condition */
bigint *bi_new()
{
    return bi_new_sized(BIGINT_MINLEAFS);
}

/* create a direct copy of bi and return it */
/* returns NULL on out of memory condition */
bigint *bi_copy(const bigint *bi)
{
    bigint *ptr;

    ptr = bi_new_sized(bi->size);
    if (ptr == NULL) return NULL;

    memcpy(ptr->data, bi->data, bi->size * sizeof(bi_leaf));
    ptr->sign = bi->sign;
    return ptr;
}

/* copy magnitude of src to previously created bigint dst */
/* returns NULL and destroys dst on out of memory condition, otherwise returns dst */
/* copying to the same bigint has no effect */
bigint *bi_copy_mag_to(bigint *dst, const bigint *src)
{
    if (dst == src)
        return dst;

    if (bi_resize(dst, src->size) == NULL) return NULL;

    memcpy(dst->data, src->data, src->size * sizeof(bi_leaf));
    memset(dst->data + src->size, 0, (dst->size - src->size) * sizeof(bi_leaf));
    return dst;
}

/* copy src to previously created bigint dst */
/* returns NULL and destroys dst on out of memory condition, otherwise returns dst */
/* copying to the same bigint has no effect */
bigint *bi_copy_to(bigint *dst, const bigint *src)
{
    if (bi_copy_mag_to(dst, src) == NULL) return NULL;
    dst->sign = src->sign;
    return dst;
}

/* assigns a previously created bigint the value 0 */
bigint *bi_clear(bigint *bi)
{
    bi->sign = 0;
    memset(bi->data, 0, bi->size * sizeof(bi_leaf));
    return bi;
}

/* assigns a previously created bigint a specified unsigned value */
bigint *bi_assignu(bigint *bi, bi_leaf value)
{
    bi_clear(bi);
    bi->data[0] = value;
    return bi;
}

/* assigns a previously created bigint a specified value */
bigint *bi_assign(bigint *bi, bi_signed_leaf value)
{
    bi_leaf uvalue;

    if (value >= 0)
        uvalue = value;
    else if (value == -value)
        uvalue = (bi_leaf) 1 << (sizeof(value)*CHAR_BIT - 1);
    else
        uvalue = (bi_leaf) -value;
    bi_assignu(bi, uvalue);
    bi->sign = value < 0;
    return bi;
}

/* assigns a previously created bigint a specified value */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_assignl(bigint *bi, bi_intmax value)
{
    bi_uintmax uvalue;

    if (value >= 0)
        uvalue = value;
    else if (-value == value)
        uvalue = (bi_uintmax) 1 << (sizeof(bi_intmax)*CHAR_BIT-1);
    else
        uvalue = -value;

    bi_assignlu(bi, uvalue);
    bi->sign = value < 0;
    return bi;
}

/* assigns a previously created bigint a specified value */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_assignlu(bigint *bi, bi_uintmax uvalue)
{
    size_t i = 0;

    bi_clear(bi);
    if (bi_resize(bi, BIGINT_LEAFS_PER_BI_INTMAX) == NULL) return NULL;

    bi->sign = 0;
    for (; i < sizeof(uvalue)*CHAR_BIT/8; i += BIGINT_LEAFBYTES)
    {
        bi->data[i/BIGINT_LEAFBYTES] = uvalue;
        uvalue >>= BIGINT_LEAFBITS;
    }

    return bi;
}

/* converts a bigint to an unsigned leaf-sized integer */
bi_leaf bi_to_intu(const bigint *bi)
{
    return bi->data[0];
}

/* converts a bigint to a signed leaf-sized integer */
bi_signed_leaf bi_to_int(const bigint *bi)
{
    bi_signed_leaf result = bi->data[0] & BIGINT_SIGNED_LEAFMAX;
    return bi->sign? -result: result;
}

/* converts a bigint to an bi_intmax */
bi_intmax bi_to_intl(const bigint *bi)
{
    bi_uintmax uvalue = bi_to_intlu(bi);

    /* TODO: may have portability issues */
    return bi->sign? -(bi_intmax) uvalue: (bi_intmax) uvalue;
}

/* converts a bigint to a bi_uintmax */
bi_uintmax bi_to_intlu(const bigint *bi)
{
    bi_uintmax uvalue = 0;
    size_t i = min(sizeof(uvalue)*CHAR_BIT/8, bi->size*BIGINT_LEAFBYTES);

    for (; i > 0; i -= BIGINT_LEAFBYTES)
    {
        uvalue <<= BIGINT_LEAFBITS;
        uvalue |= bi->data[i/BIGINT_LEAFBYTES-1];
    }

    return uvalue;
}

/* compares two bigints. returns -1 if bi < bi2, 0 if bi == b2, 1 if bi > b2 */
/* only compares magnitude, does not compare signs */
int bi_cmp_mag(const bigint *bi, const bigint *bi2)
{
    size_t i;

    for (i = bi->size; i > bi2->size; --i)
        if (bi->data[i-1] != 0)
            return 1;

    for (i = bi2->size; i > bi->size; --i)
        if (bi2->data[i-1] != 0)
            return -1;

    for (; i > 0; --i)
    {
        bi_leaf cmp = bi->data[i-1] - bi2->data[i-1];
        if (cmp)
            return cmp > bi->data[i-1]? -1: 1;
    }

    return 0;
}

/* compares bi with val. returns -1 if bi < val, 0 if bi == val, 1 if bi > val */
/* compares signs as well as magnitude */
int bi_cmp_imm(const bigint *bi, bi_signed_leaf val)
{
    size_t vsign = val < 0;
    bi_leaf v = val >= 0? (bi_leaf) val: -val == val? (bi_leaf) 1 << (sizeof(val)*CHAR_BIT - 1): (bi_leaf) -val;
    if (bi->sign > vsign) /* negative, positive */
        return -1;
    else if (bi->sign < vsign) /* positive, negative */
        return 1;
    else
    {
        int ret = bi_used(bi) > 1;

        if (!ret)
            ret = (bi->data[0] > v) - (bi->data[0] < v);

        return vsign? -ret: ret;
    }
}

/* compares bi with val. returns -1 if bi < val, 0 if bi == val, 1 if bi > val */
/* compares signs as well as magnitude */
int bi_cmp_immu(const bigint *bi, bi_leaf val)
{
    if (bi->sign) /* negative, positive */
        return -1;
    else
    {
        int ret = bi_used(bi) > 1;

        if (!ret)
            ret = (bi->data[0] > val) - (bi->data[0] < val);

        return ret;
    }
}

/* returns 1 if bi is negative, 0 if positive or zero */
int bi_is_negative(const bigint *bi)
{
    return bi->sign;
}

/* compares bi with 0. returns -1 if bi < 0, 0 if bi == 0, 1 if bi > 0 */
/* compares signs as well as magnitude */
int bi_cmp_zero(const bigint *bi)
{
    if (bi->sign)
        return -1;
    else if (bi_is_zero(bi))
        return 0;
    else
        return 1;
}

/* compares two bigints. returns -1 if bi < bi2, 0 if bi == b2, 1 if bi > b2 */
/* compares signs as well as magnitude */
int bi_cmp(const bigint *bi, const bigint *bi2)
{
    if (bi->sign > bi2->sign) /* negative, positive */
        return -1;
    else if (bi->sign < bi2->sign) /* positive, negative */
        return 1;
    else if (bi->sign == 0 && bi2->sign == 0) /* positive, positive */
        return bi_cmp_mag(bi, bi2);
    else /* negative, negative */
        return -bi_cmp_mag(bi, bi2);
}

/* returns 1 if bi is zero, zero otherwise */
int bi_is_zero(const bigint *bi)
{
    return bi_used(bi) == 0;
}

/* returns 1 if bi is positive one, zero otherwise */
int bi_is_one(const bigint *bi)
{
    return bi->data[0] == 1 && !bi->sign && bi_used(bi) == 1;
}

/* returns 1 if bi is a power of 2 (has only one bit set), zero otherwise */
int bi_is_power_of_2(const bigint *bi)
{
    size_t i = 0, cnt = 0, idx = 0;
    bi_leaf lf;

    for (; i < bi->size; ++i)
        if (bi->data[i] != 0)
        {
            if (cnt++) return 0;
            idx = i;
        }

    if (cnt != 1) return 0;

    lf = bi->data[idx];
    return (lf & (lf-1)) == 0;
}

/* returns log2(bi), rounded down */
/* returns -1 if bi is zero */
int bi_log2(const bigint *bi)
{
    return (int) bi_bitcount(bi) - 1;
}

/* returns log2(bi), rounded down */
/* returns -1 if bi is zero */
long long bi_log2l(const bigint *bi)
{
    return (long long) bi_bitcount(bi) - 1;
}

#ifdef BIGINT_ENABLE_LIBMATH
// TODO: remove floating-point dependence?
/* returns approximation of log10(bi), rounded down */
/* returns -1 if bi is zero */
/* this function will be equal to log10(bi) or log10(bi)-1 */
int bi_log10_approx(const bigint *bi)
{
    return bi_logn_approx(bi, 10);
}

/* returns approximation of log10(bi), rounded down */
/* returns -1 if bi is zero */
/* this function will be equal to log10(bi) or log10(bi)-1 */
long long bi_log10l_approx(const bigint *bi)
{
    return bi_lognl_approx(bi, 10);
}

/* returns log10(bi), rounded down */
/* returns -1 if bi is zero, -2 if out of memory */
int bi_log10(const bigint *bi)
{
    return bi_logn(bi, 10);
}

/* returns log10(bi), rounded down */
/* returns -1 if bi is zero, -2 if out of memory */
long long bi_log10l(const bigint *bi)
{
    return bi_lognl(bi, 10);
}

// TODO: remove floating-point dependence?
/* returns approximation of logn(bi), rounded down */
/* returns -1 if bi is zero or n < 2 */
/* this function will be equal to logn(bi) or logn(bi)-1 */
int bi_logn_approx(const bigint *bi, uintmax_t n)
{
    int i = bi_log2(bi);

    if (n == 2)
        return i;
    if (i < 0 || n < 2)
        return -1;
    else
        return (int) (i / log2(n));
}

/* returns approximation of logn(bi), rounded down */
/* returns -1 if bi is zero or n < 2 */
/* this function will be equal to logn(bi) or logn(bi)-1 */
long long bi_lognl_approx(const bigint *bi, uintmax_t n)
{
    long long i = bi_log2l(bi);

    if (n == 2)
        return i;
    if (i < 0 || n < 2)
        return -1;
    else
        return i / log2(n);
}

/* returns logn(bi), rounded down */
/* returns -1 if bi is zero or n < 2, -2 if out of memory */
int bi_logn(const bigint *bi, uintmax_t n)
{
    int approx = bi_logn_approx(bi, n);
    bigint *t;

    if (n < 2)
        return -1;
    else if (n == 2)
        return approx;

    t = bi_new();
    if (bi_assign(t, n) == NULL ||
        bi_uexp_assign(t, approx+1) == NULL) return -2;
    approx += bi_cmp_mag(bi, t) >= 0;
    bi_destroy(t);

    return approx;
}

/* returns logn(bi), rounded down */
/* returns -1 if bi is zero or n < 2, -2 if out of memory */
long long bi_lognl(const bigint *bi, uintmax_t n)
{
    long long approx = bi_lognl_approx(bi, n);
    bigint *t;

    if (n < 2)
        return -1;
    else if (n == 2)
        return approx;

    t = bi_new();
    if (bi_assign(t, n) == NULL ||
        bi_uexp_assign(t, approx+1) == NULL) return -2;
    approx += bi_cmp_mag(bi, t) >= 0;
    bi_destroy(t);

    return approx;
}

#endif

/* returns the number of trailing zeroes of bi */
/* returns -1 if bi is zero */
int bi_trailing_zeroes(const bigint *bi)
{
    size_t i = 0;

    for (; i < bi->size; ++i)
        if (bi->data[i] != 0)
        {
            bi_leaf l = bi->data[i];
            i *= BIGINT_LEAFBITS;
            while (!(l & 1)) ++i, l >>= 1;
            return i;
        }

    return -1;
}

/* negates bi and returns bi */
bigint *bi_negate_assign(bigint *bi)
{
    if (!bi_is_zero(bi))
        bi->sign = !bi->sign;
    else
        bi->sign = 0;

    return bi;
}

/* create a new bigint initialized to -bi and return it */
/* returns NULL on out of memory condition */
bigint *bi_negate(const bigint *bi)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_negate_assign(result);
}

/* converts bi to the absolute value of bi */
bigint *bi_abs_assign(bigint *bi)
{
    bi->sign = 0;
    return bi;
}

/* create a new bigint initialized to abs(bi) and return it */
/* returns NULL on out of memory condition */
bigint *bi_abs(const bigint *bi)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_abs_assign(result);
}

/* shift previously created bigint bi left by 1 bit */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_shl1_assign(bigint *bi)
{
    size_t i, old_size;
    bi_leaf carry = 0;

    old_size = bi->size + bi_full(bi);

    if (bi_resize(bi, old_size) == NULL)
        return NULL;

    for (i = 0; i < old_size; ++i)
    {
        bi_leaf temp = bi->data[i];
        bi->data[i] = (temp << 1) | carry;
        carry = temp >> (BIGINT_LEAFBITS - 1);
    }

    return bi;
}

/* create a new bigint initialized to bi << 1 and return it */
/* returns NULL on out of memory condition */
bigint *bi_shl1(const bigint *bi)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_shl1_assign(result);
}

/* shift previously created bigint bi left by n leaves */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_shl_leaf_assign(bigint *bi, size_t n)
{
    size_t old_size = bi_used(bi);

    if (bi_resize(bi, old_size + n) == NULL)
        return NULL;

    if (n)
    {
        memmove(bi->data + n, bi->data, old_size * sizeof(bi_leaf));
        memset(bi->data, 0, n * sizeof(bi_leaf));
    }

    return bi;
}

/* create a new bigint initialized to bi << n (leaves) and return it */
/* returns NULL on out of memory condition */
bigint *bi_shl_leaf(const bigint *bi, size_t n)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_shl_leaf_assign(result, n);
}

/* shift previously created bigint bi right by n leaves, rounding toward 0 */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_shr_leaf_assign(bigint *bi, size_t n)
{
    size_t t = min(n, bi->size);

    if (t)
    {
        memmove(bi->data, bi->data + t, (bi->size - t) * sizeof(bi_leaf));
        memset(bi->data + (bi->size - t), 0, t * sizeof(bi_leaf));
    }

    return bi;
}

/* create a new bigint initialized to bi >> n (leaves) and return it, rounding toward 0 */
/* returns NULL on out of memory condition */
bigint *bi_shr_leaf(const bigint *bi, size_t n)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_shr_leaf_assign(result, n);
}

/* shift previously created bigint bi left by n bits */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_shl_assign(bigint *bi, size_t n)
{
    size_t t, i, old_size;
    bi_leaf carry = 0;

    old_size = bi_used(bi) + 1;
    t = n / BIGINT_LEAFBITS;

    if (bi_resize(bi, old_size + t) == NULL)
        return NULL;

    if (t)
    {
        memmove(bi->data + t, bi->data, old_size * sizeof(bi_leaf));
        memset(bi->data, 0, t * sizeof(bi_leaf));
    }

    t = n % BIGINT_LEAFBITS;
    if (t)
    {
        for (i = 0; i < bi->size; ++i)
        {
            bi_leaf temp = bi->data[i];
            bi->data[i] = (temp << t) | carry;
            carry = temp >> (BIGINT_LEAFBITS - t);
        }
    }

    return bi;
}

/* create a new bigint initialized to bi << n and return it */
/* returns NULL on out of memory condition */
bigint *bi_shl(const bigint *bi, size_t n)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_shl_assign(result, n);
}

/* shift previously created bigint bi right by n bits, rounding toward 0 */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_shr_assign(bigint *bi, size_t n)
{
    size_t t, i;
    bi_leaf carry = 0;

    t = min(n / BIGINT_LEAFBITS, bi->size);

    if (t)
    {
        memmove(bi->data, bi->data + t, (bi->size - t) * sizeof(bi_leaf));
        memset(bi->data + (bi->size - t), 0, t * sizeof(bi_leaf));
    }

    t = n % BIGINT_LEAFBITS;
    if (t)
    {
        for (i = bi->size; i > 0; --i)
        {
            bi_leaf temp = bi->data[i-1];
            bi->data[i-1] = (temp >> t) | carry;
            carry = temp << (BIGINT_LEAFBITS - t);
        }
    }

    return bi;
}

/* create a new bigint initialized to bi >> n and return it, rounding toward 0 */
/* returns NULL on out of memory condition */
bigint *bi_shr(const bigint *bi, size_t n)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_shr_assign(result, n);
}

/* create a new bigint initialized to bi & bi2 and return it */
/* returns NULL on out of memory condition */
bigint *bi_and(const bigint *bi, const bigint *bi2)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_and_assign(result, bi2);
}

/* bitwise ANDs previously created bigint bi with bi2; the sign is not affected */
bigint *bi_and_assign(bigint *bi, const bigint *bi2)
{
    size_t i = 0, m = min(bi->size, bi2->size);
    memset(bi->data + m, 0, (bi->size - m) * sizeof(bi_leaf));
    for (; i < m; ++i)
        bi->data[i] &= bi2->data[i];
    return bi;
}

/* create a new bigint initialized to bi | bi2 and return it */
/* returns NULL on out of memory condition */
bigint *bi_or(const bigint *bi, const bigint *bi2)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_or_assign(result, bi2);
}

/* bitwise ORs previously created bigint bi with bi2; the sign is not affected */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_or_assign(bigint *bi, const bigint *bi2)
{
    size_t i = 0, m = max(bi->size, bi2->size), bi2_size = bi2->size;
    if (bi_resize(bi, m) == NULL) return NULL;
    for (; i < bi2_size; ++i)
        bi->data[i] |= bi2->data[i];
    return bi;
}

/* create a new bigint initialized to bi ^ bi2 and return it */
/* returns NULL on out of memory condition */
bigint *bi_xor(const bigint *bi, const bigint *bi2)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_xor_assign(result, bi2);
}

/* bitwise XORs previously created bigint bi with bi2; the sign is not affected */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_xor_assign(bigint *bi, const bigint *bi2)
{
    size_t i = 0, m = max(bi->size, bi2->size), bi2_size = bi2->size;
    if (bi_resize(bi, m) == NULL) return NULL;
    for (; i < bi2_size; ++i)
        bi->data[i] ^= bi2->data[i];
    return bi;
}

static bigint *bi_add_invert(bigint *dst, const bigint *bi, const bigint *bi2, int invert2);
static bigint *bi_add_imm_invert(bigint *dst, const bigint *bi, bi_signed_leaf val, int invert2);
static bigint *bi_sub_invert(bigint *dst, const bigint *bi, const bigint *bi2, int invert2);
static bigint *bi_sub_imm_invert(bigint *dst, const bigint *bi, bi_signed_leaf val, int invert2);

/* adds two bigints together and returns the result in dst, or a new bigint if dst is NULL */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_add_invert(bigint *dst, const bigint *bi, const bigint *bi2, int invert2)
{
    bigint *result = dst;
    size_t i = 0, bi_size = bi->size - !bi_full(bi), bi2_size = bi2->size - !bi_full(bi2), end;
    bi_dleaf carry = 0;
    int invert1;

    invert1 = bi->sign;
    invert2 = invert2? !bi2->sign: bi2->sign;

    /* add only if signs are equal, subtract otherwise */
    if (invert1 != invert2)
    {
        /* negative + positive = negative - negative */
        /* positive + negative = positive - positive */
        return bi_sub_invert(dst, bi, bi2, 1);
    }

    if (dst == NULL)
    {
        result = bi_new_sized(max(BIGINT_MINLEAFS, max(bi_size, bi2_size) + 1));
        if (result == NULL) return NULL;
    }
    else
    {
        if (bi_resize(result, max(bi_size, bi2_size) + 1) == NULL)
            return NULL;
    }

    end = min(bi_size, bi2_size);
    for (; i < end; ++i)
    {
        carry += bi->data[i];
        carry += bi2->data[i];
        result->data[i] = carry;
        carry >>= BIGINT_LEAFBITS;
    }

    if (bi_size > bi2_size)
        bi2 = bi, bi2_size = bi_size;

    for (; i < bi2_size; ++i) /* add the rest of the larger number */
    {
        carry += bi2->data[i];
        result->data[i] = carry;
        carry >>= BIGINT_LEAFBITS;
    }

    result->data[i] = carry;
    result->sign = invert1;

    return result;
}

/* adds val to bi and returns the result in dst, or a new bigint if dst is NULL */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_add_imm_invert(bigint *dst, const bigint *bi, bi_signed_leaf val, int invert2)
{
    bigint *result = dst;
    size_t i = 0, bi_size = bi->size - !bi_full(bi);
    bi_dleaf carry = val >= 0? (bi_dleaf) val: -val == val? (bi_dleaf) 1 << (sizeof(val)*CHAR_BIT - 1): (bi_dleaf) -val;
    int invert1;

    invert1 = bi->sign;
    invert2 = invert2? val >= 0: val < 0;

    /* add only if signs are equal, subtract otherwise */
    if (invert1 != invert2)
    {
        /* negative + positive = negative - negative */
        /* positive + negative = positive - positive */
        return bi_sub_imm_invert(dst, bi, val, 1);
    }

    if (dst == NULL)
    {
        result = bi_new_sized(max(BIGINT_MINLEAFS, bi_size + 1));
        if (result == NULL) return NULL;
    }
    else
    {
        if (bi_resize(result, bi_size + 1) == NULL)
            return NULL;
    }

    for (; i < bi_size; ++i)
    {
        carry += bi->data[i];
        result->data[i] = carry;
        carry >>= BIGINT_LEAFBITS;
    }

    result->data[i] = carry;
    result->sign = invert1;

    return result;
}

/* subtracts bi2 from bi and returns the result in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_sub_invert(bigint *dst, const bigint *bi, const bigint *bi2, int invert2)
{
    bigint *result = dst;
    size_t i = 0, bi_size, bi2_size, end;
    bi_dleaf borrow = 0;
    int invert1;
    int mag_cmp = 0, swapped = 0;

    invert1 = bi->sign;
    invert2 = invert2? !bi2->sign: bi2->sign;

    mag_cmp = bi_cmp_mag(bi, bi2);

    if (mag_cmp < 0)
    {
        const bigint *temp = bi;
        bi = bi2;
        bi2 = temp;

        swapped = invert1;
        invert1 = invert2;
        invert2 = swapped;

        swapped = 1;
    }
    bi_size = bi->size;
    bi2_size = bi2->size;

    /* subtract only if signs are equal, add otherwise */
    if (invert1 != invert2)
    {
        /* negative - positive = negative + negative */
        /* positive - negative = positive + positive */
        return bi_add_invert(dst, bi, bi2, 1);
    }

    if (dst == NULL)
    {
        result = bi_new_sized(max(BIGINT_MINLEAFS, bi_size));
        if (result == NULL) return NULL;
    }
    else
    {
        if (bi_resize(result, bi_size) == NULL)
            return NULL;
    }

    end = min(bi_size, bi2_size);
    for (; i < end; ++i)
    {
        borrow = bi->data[i] - borrow;
        borrow -= bi2->data[i];
        result->data[i] = borrow;
        borrow >>= BIGINT_LEAFBITS*2-1;
    }

    for (; i < bi_size; ++i) /* subtract the rest (bi is guaranteed to be larger than bi2) */
    {
        borrow = bi->data[i] - borrow;
        result->data[i] = borrow;
        borrow >>= BIGINT_LEAFBITS*2-1;
    }

    result->sign = (invert1? !swapped: swapped) && mag_cmp;

    return result;
}

/* subtracts bi2 from bi, and sets bi to the result */
/* assumes bi2 is less than or equal to bi */
static bigint *bi_fast_sub_assign(bigint *bi, const bigint *bi2)
{
    size_t i = 0, end = bi->size < bi2->size? bi->size: bi2->size;
    bi_dleaf borrow = 0;

    for (; i < end; ++i)
    {
        borrow = bi->data[i] - borrow;
        borrow -= bi2->data[i];
        bi->data[i] = borrow;
        borrow >>= BIGINT_LEAFBITS*2-1;
    }

    for (; i < bi->size; ++i) /* subtract the rest */
    {
        borrow = bi->data[i] - borrow;
        bi->data[i] = borrow;
        borrow >>= BIGINT_LEAFBITS*2-1;
    }

    return bi;
}

/* subtracts bi2 from bi and returns the result in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_sub_imm_invert(bigint *dst, const bigint *bi, bi_signed_leaf val, int invert2)
{
    bigint *result = dst;
    size_t i = 0, bi_size;
    bi_dleaf borrow = val >= 0? (bi_dleaf) val: -val == val? (bi_dleaf) 1 << (sizeof(val)*CHAR_BIT - 1): (bi_dleaf) -val;
    int invert1;
    int mag_cmp = 0, swapped = 0;

    invert1 = bi->sign;
    invert2 = invert2? val >= 0: val < 0;

    mag_cmp = (bi->data[0] > borrow || bi_used(bi) > 1) - (bi->data[0] < borrow);

    if (mag_cmp < 0)
    {
        bi_leaf temp = bi->data[0];
        bi->data[0] = borrow;
        borrow = temp;

        swapped = invert1;
        invert1 = invert2;
        invert2 = swapped;

        swapped = 1;
    }
    bi_size = bi->size - !bi_full(bi);

    /* subtract only if signs are equal, add otherwise */
    if (invert1 != invert2)
    {
        /* negative - positive = negative + negative */
        /* positive - negative = positive + positive */
        return bi_add_imm_invert(dst, bi, val, 1);
    }

    if (dst == NULL)
    {
        result = bi_new_sized(max(BIGINT_MINLEAFS, bi_size));
        if (result == NULL) return NULL;
    }

    for (; i < result->size; ++i)
    {
        borrow = bi->data[i] - borrow;
        result->data[i] = borrow;
        borrow >>= BIGINT_LEAFBITS*2-1;
    }

    result->sign = (invert1? !swapped: swapped) && mag_cmp;

    return result;
}

/* adds two bigints together and returns the result in a new bigint */
/* returns NULL on out of memory or erroneous parameter conditions; the parameters are not destroyed */
bigint *bi_add(const bigint *bi, const bigint *bi2)
{
    return bi_add_invert(NULL, bi, bi2, 0);
}

/* adds bi and bi2 and sets bi to the result. returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_add_assign(bigint *bi, const bigint *bi2)
{
    return bi_add_invert(bi, bi, bi2, 0);
}

/* adds bi and val and returns the result in a new bigint */
/* returns NULL on out of memory or erroneous parameter conditions; the parameters are not destroyed */
bigint *bi_add_immediate(const bigint *bi, bi_signed_leaf val)
{
    return bi_add_imm_invert(NULL, bi, val, 0);
}

/* adds bi and bi2 and sets bi to the result. returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_add_immediate_assign(bigint *bi, bi_signed_leaf val)
{
    return bi_add_imm_invert(bi, bi, val, 0);
}

/* subtracts bi2 from bi and returns the result in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
bigint *bi_sub(const bigint *bi, const bigint *bi2)
{
    return bi_sub_invert(NULL, bi, bi2, 0);
}

/* subtracts bi2 from bi and sets bi to the result. returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_sub_assign(bigint *bi, const bigint *bi2)
{
    return bi_sub_invert(bi, bi, bi2, 0);
}

/* subtracts val from bi and returns the result in a new bigint */
/* returns NULL on out of memory or erroneous parameter conditions; the parameters are not destroyed */
bigint *bi_sub_immediate(const bigint *bi, bi_signed_leaf val)
{
    return bi_sub_imm_invert(NULL, bi, val, 0);
}

/* subtracts val from bi and sets bi to the result. returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_sub_immediate_assign(bigint *bi, bi_signed_leaf val)
{
    return bi_sub_imm_invert(bi, bi, val, 0);
}

/* multiplies bi and bi2 using standard multiplication */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_mul_internal(const bigint *bi, const bigint *bi2, size_t sz, size_t sz2)
{
    bi_dleaf low_carry = 0, high_carry = 0;
    bigint *result;

    size_t i, j;

    if (sz == 0 || sz2 == 0) return bi_new();
    else if (sz == 1) {return bi_mul_immediateu(bi2, bi->data[0]);}
    else if (sz2 == 1) {return bi_mul_immediateu(bi, bi2->data[0]);}

    result = bi_new_sized(max(sz + sz2, BIGINT_MINLEAFS));
    if (result == NULL) return NULL;

    for (i = 0; i < sz + sz2; ++i)
    {
        size_t j_start, j_end;

        j_start = i >= sz2? i - sz2 + 1: 0;
        j_end = min(i, sz-1);

        for (j = j_start; j <= j_end; ++j)
        {
            bi_dleaf a, b, save;

            a = bi->data[j];
            b = bi2->data[i - j];

            save = low_carry;
            low_carry += a * b;
            high_carry += low_carry < save;
        }

        result->data[i] = low_carry;
        low_carry = (low_carry >> BIGINT_LEAFBITS) | (high_carry << BIGINT_LEAFBITS);
        high_carry >>= BIGINT_LEAFBITS;
    }

    return result;
}

#ifdef BIGINT_ENABLE_PTHREADS
static bigint *bi_karatsuba_mul(const bigint *bi, const bigint *bi2, size_t bi_used_leafs, size_t bi2_used_leafs, int use_pthreads);
static bigint *bi_toom_cook_mul(const bigint *bi, const bigint *bi2, size_t bi_used_leafs, size_t bi2_used_leafs, int use_pthreads);

struct bigint_pthread_metadata
{
    bigint *op1, *op2;
    size_t used1, used2;
};

void *bi_karatsuba_pthread(void *p)
{
    struct bigint_pthread_metadata *mdata = (struct bigint_pthread_metadata *) p;
    return bi_karatsuba_mul(mdata->op1, mdata->op2, mdata->used1, mdata->used2, 0);
}

void *bi_toom_cook_pthread(void *p)
{
    struct bigint_pthread_metadata *mdata = (struct bigint_pthread_metadata *) p;
    return bi_toom_cook_mul(mdata->op1, mdata->op2, mdata->used1, mdata->used2, 0);
}
#endif

/* multiplies bi and bi2 using Karatsuba multiplication, reverting to standard multiplication with small enough chunks */
/* returns NULL on out of memory condition; the parameters are not destroyed */
#ifdef BIGINT_ENABLE_PTHREADS
static bigint *bi_karatsuba_mul(const bigint *bi, const bigint *bi2, size_t bi_used_leafs, size_t bi2_used_leafs, int use_pthreads)
#else
static bigint *bi_karatsuba_mul(const bigint *bi, const bigint *bi2, size_t bi_used_leafs, size_t bi2_used_leafs)
#endif
{
    if (bi_used_leafs < BIGINT_KARATSUBA_MIN_LEAFS || bi2_used_leafs < BIGINT_KARATSUBA_MIN_LEAFS)
        return bi_mul_internal(bi, bi2, bi_used_leafs, bi2_used_leafs);
    else
    {
        size_t m = max(bi_used_leafs, bi2_used_leafs) / 2 + 1;
        size_t used1, used2;
        bigint high1, low1, high2, low2, *z0 = NULL, *z1 = NULL, *z2 = NULL, *temp = NULL, *temp2 = NULL;

        low1 = high1 = *bi;
        low1.sign = high1.sign = 0;
        low1.size = min(m, bi_used_leafs);
        if (bi_used_leafs > m)
        {
            high1.data += m;
            high1.size = bi_used_leafs - m;
        }
        else
            high1.size = 0;
        low2 = high2 = *bi2;
        low2.sign = high2.sign = 0;
        low2.size = min(m, bi2_used_leafs);
        if (bi2_used_leafs > m)
        {
            high2.data += m;
            high2.size = bi2_used_leafs - m;
        }
        else
            high2.size = 0;

#ifdef BIGINT_ENABLE_PTHREADS
        if (use_pthreads)
        {
            void *vz0, *vz1, *vz2;
            pthread_t low_thread, mid_thread, high_thread;
            struct bigint_pthread_metadata low_data, mid_data, high_data;

            low_data.op1 = &low1;
            low_data.op2 = &low2;
            low_data.used1 = bi_used(&low1);
            low_data.used2 = bi_used(&low2);

            if (high1.size)
                temp = bi_add(&low1, &high1);
            else
                temp = &low1;
            if (temp == NULL) goto cleanup;
            if (high2.size)
                temp2 = bi_add(&low2, &high2);
            else
                temp2 = &low2;
            if (temp2 == NULL) goto cleanup;

            mid_data.op1 = temp;
            mid_data.op2 = temp2;
            mid_data.used1 = bi_used(temp);
            mid_data.used2 = bi_used(temp2);

            high_data.op1 = &high1;
            high_data.op2 = &high2;
            high_data.used1 = bi_used(&high1);
            high_data.used2 = bi_used(&high2);

            if (pthread_create(&low_thread, NULL, bi_karatsuba_pthread, &low_data)) goto cleanup;
            if (pthread_create(&mid_thread, NULL, bi_karatsuba_pthread, &mid_data)) goto cleanup;
            if (pthread_create(&high_thread, NULL, bi_karatsuba_pthread, &high_data)) goto cleanup;

            if (pthread_join(low_thread, &vz0)) goto cleanup;
            z0 = vz0;
            if (z0 == NULL) goto cleanup;

            if (pthread_join(mid_thread, &vz1)) goto cleanup;
            z1 = vz1;
            if (z1 == NULL) goto cleanup;

            if (pthread_join(high_thread, &vz2)) goto cleanup;
            z2 = vz2;
            if (z2 == NULL) goto cleanup;

            if (temp != &low1) bi_destroy(temp); temp=NULL;
            if (temp2 != &low2) bi_destroy(temp2); temp2=NULL;
        }
        else
        {
            used1 = bi_used(&low1);
            used2 = bi_used(&low2);
            if (used1 && used2)
            {
                z0 = bi_karatsuba_mul(&low1, &low2, used1, used2, 0);
                if (z0 == NULL) goto cleanup;
            }
            /* calculate z1 = (low1 + high1) * (low2 + high2) */
            if (high1.size)
                temp = bi_add(&low1, &high1);
            else
                temp = &low1;
            if (temp == NULL) goto cleanup;
            if (high2.size)
                temp2 = bi_add(&low2, &high2);
            else
                temp2 = &low2;
            if (temp2 == NULL) goto cleanup;
            z1 = bi_karatsuba_mul(temp, temp2, bi_used(temp), bi_used(temp2), 0);
            if (z1 == NULL) goto cleanup;
            if (temp != &low1) bi_destroy(temp); temp=NULL;
            if (temp2 != &low2) bi_destroy(temp2); temp2=NULL;
            /* calculate z2 = high1 * high2 */
            z2 = bi_karatsuba_mul(&high1, &high2, bi_used(&high1), bi_used(&high2), 0);
            if (z2 == NULL) goto cleanup;
        }
        /* calculate z1 = z1 - z2 - z0 */
        bi_fast_sub_assign(z1, z2);
        if (z0 != NULL) bi_fast_sub_assign(z1, z0);
        /* calculate z2 = z2 * 2^(2*m) + z1 * 2^m + z0 */
        if ((z2 = bi_shl_leaf_assign(z2, m)) == NULL) goto cleanup;
        if ((z2 = bi_add_assign(z2, z1)) == NULL) goto cleanup;
        bi_destroy(z1); z1=NULL;
        if ((z2 = bi_shl_leaf_assign(z2, m)) == NULL) goto cleanup;
        if (z0 != NULL && (z2 = bi_add_assign(z2, z0)) == NULL) goto cleanup;
        bi_destroy(z0);
#else
        used1 = bi_used(&low1);
        used2 = bi_used(&low2);
        if (used1 && used2)
        {
            z0 = bi_karatsuba_mul(&low1, &low2, used1, used2);
            if (z0 == NULL) goto cleanup;
        }
        /* calculate z1 = (low1 + high1) * (low2 + high2) */
        if (high1.size)
            temp = bi_add(&low1, &high1);
        else
            temp = &low1;
        if (temp == NULL) goto cleanup;
        if (high2.size)
            temp2 = bi_add(&low2, &high2);
        else
            temp2 = &low2;
        if (temp2 == NULL) goto cleanup;
        z1 = bi_karatsuba_mul(temp, temp2, bi_used(temp), bi_used(temp2));
        if (z1 == NULL) goto cleanup;
        if (temp != &low1) bi_destroy(temp); temp=NULL;
        if (temp2 != &low2) bi_destroy(temp2); temp2=NULL;
        /* calculate z2 = high1 * high2 */
        z2 = bi_karatsuba_mul(&high1, &high2, bi_used(&high1), bi_used(&high2));
        if (z2 == NULL) goto cleanup;
        /* calculate z1 = z1 - z2 - z0 */
        bi_fast_sub_assign(z1, z2);
        if (z0 != NULL) bi_fast_sub_assign(z1, z0);
        /* calculate z2 = z2 * 2^(2*m) + z1 * 2^m + z0 */
        if ((z2 = bi_shl_leaf_assign(z2, m)) == NULL) goto cleanup;
        if ((z2 = bi_add_assign(z2, z1)) == NULL) goto cleanup;
        bi_destroy(z1); z1=NULL;
        if ((z2 = bi_shl_leaf_assign(z2, m)) == NULL) goto cleanup;
        if (z0 != NULL && (z2 = bi_add_assign(z2, z0)) == NULL) goto cleanup;
        bi_destroy(z0);
#endif

        z2->sign = bi->sign != bi2->sign;
        return z2;
cleanup:
        bi_destroy(z0);
        if (temp != &low1) bi_destroy(temp);
        if (temp2 != &low2) bi_destroy(temp2);
        bi_destroy(z1);
        bi_destroy(z2);
        return NULL;
    }
}

/* multiplies bi and bi2 using Toom-3 multiplication, reverting to standard multiplication with small enough chunks */
/* returns NULL on out of memory condition; the parameters are not destroyed */
#ifdef BIGINT_ENABLE_PTHREADS
static bigint *bi_toom_cook_mul(const bigint *bi, const bigint *bi2, size_t bi_used_leafs, size_t bi2_used_leafs, int use_pthreads)
{
    if (bi_used_leafs < BIGINT_TOOM_COOK_MIN_LEAFS || bi2_used_leafs < BIGINT_TOOM_COOK_MIN_LEAFS)
        return bi_karatsuba_mul(bi, bi2, bi_used_leafs, bi2_used_leafs, 0);
    else
    {
        size_t m = max(bi_used_leafs, bi2_used_leafs) / 3 + 1;
        size_t two_m = 2 * m;
        bigint high1, mid1, low1, high2, mid2, low2;
        bigint *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL, *temp = NULL, *temp2 = NULL;
        bigint *temp3 = NULL, *temp4 = NULL, *temp5 = NULL, *temp6 = NULL, *intermediate1 = NULL, *intermediate2 = NULL;
        bigint *r0 = NULL, *r1 = NULL, *r2 = NULL, *r3 = NULL, *r4 = NULL;

        low1 = mid1 = high1 = *bi;
        low1.sign = mid1.sign = high1.sign = 0;
        low1.size = min(m, bi_used_leafs);
        if (bi_used_leafs > m)
        {
            mid1.data += m;
            mid1.size = min(bi_used_leafs - m, m);
            if (bi_used_leafs > two_m)
            {
                high1.data += two_m;
                high1.size = bi_used_leafs - two_m;
            }
            else
                high1.size = 0;
        }
        else
            mid1.size = high1.size = 0;

        low2 = mid2 = high2 = *bi2;
        low2.sign = mid2.sign = high2.sign = 0;
        low2.size = min(m, bi2_used_leafs);
        if (bi2_used_leafs > m)
        {
            mid2.data += m;
            mid2.size = min(bi2_used_leafs - m, m);
            if (bi2_used_leafs > two_m)
            {
                high2.data += two_m;
                high2.size = bi2_used_leafs - two_m;
            }
            else
                high2.size = 0;
        }
        else
            mid2.size = high2.size = 0;

        /* for each quadratic polynomial (hx^2 + mx + l), using points (0, 1, -1, -2, infinity):
         *   p0 = l
         *   p1 = h + m + l
         *   p2 = h - m + l
         *   p3 = 4h - 2m + l
         *   p4 = h
         */

        /* calculate intermediate1 = h + l */
        if (high1.size)
            intermediate1 = bi_add(&high1, &low1);
        else
            intermediate1 = &low1;

        /* calculate intermediate2 = h + l */
        if (high2.size)
            intermediate2 = bi_add(&high2, &low2);
        else
            intermediate2 = &low2;

        if (use_pthreads)
        {
            void *vp0, *vp1, *vp2, *vp3, *vp4;
            pthread_t p0_thread, p1_thread, p2_thread, p3_thread, p4_thread;
            struct bigint_pthread_metadata p0_data, p1_data, p2_data, p3_data, p4_data;

            // P0
            p0_data.op1 = &low1;
            p0_data.op2 = &low2;
            p0_data.used1 = bi_used(&low1);
            p0_data.used2 = bi_used(&low2);

            // P1
            temp = bi_add(intermediate1, &mid1);
            if (temp == NULL) goto cleanup;

            temp2 = bi_add(intermediate2, &mid2);
            if (temp2 == NULL) goto cleanup;

            p1_data.op1 = temp;
            p1_data.op2 = temp2;
            p1_data.used1 = bi_used(temp);
            p1_data.used2 = bi_used(temp2);

            // P2
            temp3 = bi_sub(intermediate1, &mid1);
            if (temp3 == NULL) goto cleanup;

            temp4 = bi_sub(intermediate2, &mid2);
            if (temp4 == NULL) goto cleanup;

            p2_data.op1 = temp3;
            p2_data.op2 = temp4;
            p2_data.used1 = bi_used(temp3);
            p2_data.used2 = bi_used(temp4);

            // P3
            temp5 = bi_copy(temp3);
            if (temp5 == NULL) goto cleanup;

            temp6 = bi_copy(temp4);
            if (temp6 == NULL) goto cleanup;

            if ((temp5 = bi_add_assign(temp5, &high1)) == NULL) goto cleanup;
            if ((temp5 = bi_shl1_assign(temp5)) == NULL) goto cleanup;
            if ((temp5 = bi_sub_assign(temp5, &low1)) == NULL) goto cleanup;

            if ((temp6 = bi_add_assign(temp6, &high2)) == NULL) goto cleanup;
            if ((temp6 = bi_shl1_assign(temp6)) == NULL) goto cleanup;
            if ((temp6 = bi_sub_assign(temp6, &low2)) == NULL) goto cleanup;

            p3_data.op1 = temp5;
            p3_data.op2 = temp6;
            p3_data.used1 = bi_used(temp5);
            p3_data.used2 = bi_used(temp6);

            // P4
            p4_data.op1 = &high1;
            p4_data.op2 = &high2;
            p4_data.used1 = bi_used(&high1);
            p4_data.used2 = bi_used(&high2);

            if (pthread_create(&p0_thread, NULL, bi_toom_cook_pthread, &p0_data)) goto cleanup;
            if (pthread_create(&p1_thread, NULL, bi_toom_cook_pthread, &p1_data)) goto cleanup;
            if (pthread_create(&p2_thread, NULL, bi_toom_cook_pthread, &p2_data)) goto cleanup;
            if (pthread_create(&p3_thread, NULL, bi_toom_cook_pthread, &p3_data)) goto cleanup;
            if (pthread_create(&p4_thread, NULL, bi_toom_cook_pthread, &p4_data)) goto cleanup;

            if (pthread_join(p0_thread, &vp0)) goto cleanup;
            p0 = vp0;
            if (p0 == NULL) goto cleanup;

            if (pthread_join(p1_thread, &vp1)) goto cleanup;
            p1 = vp1;
            if (p1 == NULL) goto cleanup;

            if (pthread_join(p2_thread, &vp2)) goto cleanup;
            p2 = vp2;
            if (p2 == NULL) goto cleanup;

            if (pthread_join(p3_thread, &vp3)) goto cleanup;
            p3 = vp3;
            if (p3 == NULL) goto cleanup;

            if (pthread_join(p4_thread, &vp4)) goto cleanup;
            p4 = vp4;
            if (p4 == NULL) goto cleanup;

            bi_destroy(temp); temp=NULL;
            bi_destroy(temp2); temp2=NULL;
            bi_destroy(temp3); temp3=NULL;
            bi_destroy(temp4); temp4=NULL;
            bi_destroy(temp5); temp5=NULL;
            bi_destroy(temp6); temp6=NULL;
            if (intermediate1 != &low1) bi_destroy(intermediate1); intermediate1=NULL;
            if (intermediate2 != &low2) bi_destroy(intermediate2); intermediate2=NULL;
        }
        else
        {
            /* calculate p0 = low1 * low2 */
            p0 = bi_toom_cook_mul(&low1, &low2, bi_used(&low1), bi_used(&low2), 0);
            if (p0 == NULL) goto cleanup;

            /* calculate p1 = (low1 + mid1 + high1) * (low2 + mid2 + high2) */
            temp = bi_add(intermediate1, &mid1);
            if (temp == NULL) goto cleanup;

            temp2 = bi_add(intermediate2, &mid2);
            if (temp2 == NULL) goto cleanup;

            p1 = bi_toom_cook_mul(temp, temp2, bi_used(temp), bi_used(temp2), 0);
            if (p1 == NULL) goto cleanup;
            bi_destroy(temp); temp=NULL;
            bi_destroy(temp2); temp2=NULL;

            /* calculate p2 = (low1 - mid1 + high1) * (low2 - mid2 + high2) */
            temp = bi_sub(intermediate1, &mid1);
            if (temp == NULL) goto cleanup;

            temp2 = bi_sub(intermediate2, &mid2);
            if (temp2 == NULL) goto cleanup;

            p2 = bi_toom_cook_mul(temp, temp2, bi_used(temp), bi_used(temp2), 0);
            if (p2 == NULL) goto cleanup;

            if (intermediate1 != &low1) bi_destroy(intermediate1); intermediate1=NULL;
            if (intermediate2 != &low2) bi_destroy(intermediate2); intermediate2=NULL;
            // Do not cleanup temp and temp2 because the next step needs them

            /* calculate p3 = (low1 - 2*mid1 + 4*high1) * (low2 - 2*mid2 + 4*high2) */
            if ((temp = bi_add_assign(temp, &high1)) == NULL) goto cleanup;
            if ((temp = bi_shl1_assign(temp)) == NULL) goto cleanup;
            if ((temp = bi_sub_assign(temp, &low1)) == NULL) goto cleanup;

            if ((temp2 = bi_add_assign(temp2, &high2)) == NULL) goto cleanup;
            if ((temp2 = bi_shl1_assign(temp2)) == NULL) goto cleanup;
            if ((temp2 = bi_sub_assign(temp2, &low2)) == NULL) goto cleanup;

            p3 = bi_toom_cook_mul(temp, temp2, bi_used(temp), bi_used(temp2), 0);
            if (p3 == NULL) goto cleanup;
            bi_destroy(temp); temp=NULL;
            bi_destroy(temp2); temp2=NULL;

            /* calculate p4 = high1 * high2 */
            p4 = bi_toom_cook_mul(&high1, &high2, bi_used(&high1), bi_used(&high2), 0);
            if (p4 == NULL) goto cleanup;
        }

        /* calculate interpolation */
        r0 = p0;
        r4 = p4;

        r3 = bi_sub(p3, p1);
        bi_destroy(p3); p3=NULL; // p3 no longer needed
        if (r3 == NULL) goto cleanup;
        if ((r3 = bi_div_3_assign(r3)) == NULL) goto cleanup;

        r1 = bi_sub(p1, p2);
        if (r1 == NULL) goto cleanup;
        (void) bi_shr_assign(r1, 1);
        bi_destroy(p1); p1=NULL; // p1 no longer needed

        r2 = bi_sub(p2, p0);
        if (r2 == NULL) goto cleanup;
        bi_destroy(p2); p2=NULL; // p2 no longer needed

        (void) bi_negate_assign(r3);
        if ((r3 = bi_add_assign(r3, r2)) == NULL) goto cleanup;
        (void) bi_shr_assign(r3, 1);

        temp = bi_shl1(r4);
        if ((r3 = bi_add_assign(r3, temp)) == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;

        if ((r2 = bi_add_assign(r2, r1)) == NULL) goto cleanup;
        if ((r2 = bi_sub_assign(r2, r4)) == NULL) goto cleanup;

        if ((r1 = bi_sub_assign(r1, r3)) == NULL) goto cleanup;

        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r3)) == NULL) goto cleanup;
        bi_destroy(r3); r3=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r2)) == NULL) goto cleanup;
        bi_destroy(r2); r2=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r1)) == NULL) goto cleanup;
        bi_destroy(r1); r1=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r0)) == NULL) goto cleanup;
        bi_destroy(r0); r0=NULL;

        // no need to free p0 because it's an alias for r0
        // no need to free p4 because we are returning it as r4 (r4 is an alias)
        // p1, p2, p3, temp, temp2, intermediate1, and intermediate2 are already freed
        r4->sign = bi->sign != bi2->sign;
        return r4;

cleanup:
        bi_destroy(r1);
        bi_destroy(r2);
        bi_destroy(r3);

        bi_destroy(p0);
        bi_destroy(temp);
        bi_destroy(temp2);
        bi_destroy(temp3);
        bi_destroy(temp4);
        bi_destroy(temp5);
        bi_destroy(temp6);
        if (intermediate1 != &low1) bi_destroy(intermediate1);
        if (intermediate2 != &low2) bi_destroy(intermediate2);
        bi_destroy(p1);
        bi_destroy(p2);
        bi_destroy(p3);
        bi_destroy(p4);
        return NULL;
    }
}
#else
static bigint *bi_toom_cook_mul(const bigint *bi, const bigint *bi2, size_t bi_used_leafs, size_t bi2_used_leafs)
{
    if (bi_used_leafs < BIGINT_TOOM_COOK_MIN_LEAFS || bi2_used_leafs < BIGINT_TOOM_COOK_MIN_LEAFS)
        return bi_karatsuba_mul(bi, bi2, bi_used_leafs, bi2_used_leafs);
    else
    {
        size_t m = max(bi_used_leafs, bi2_used_leafs) / 3 + 1;
        size_t two_m = 2 * m;
        bigint high1, mid1, low1, high2, mid2, low2;
        bigint *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL, *temp = NULL, *temp2 = NULL, *intermediate1 = NULL, *intermediate2 = NULL;
        bigint *r0 = NULL, *r1 = NULL, *r2 = NULL, *r3 = NULL, *r4 = NULL;

        low1 = mid1 = high1 = *bi;
        low1.sign = mid1.sign = high1.sign = 0;
        low1.size = min(m, bi_used_leafs);
        if (bi_used_leafs > m)
        {
            mid1.data += m;
            mid1.size = min(bi_used_leafs - m, m);
            if (bi_used_leafs > two_m)
            {
                high1.data += two_m;
                high1.size = bi_used_leafs - two_m;
            }
            else
                high1.size = 0;
        }
        else
            mid1.size = high1.size = 0;

        low2 = mid2 = high2 = *bi2;
        low2.sign = mid2.sign = high2.sign = 0;
        low2.size = min(m, bi2_used_leafs);
        if (bi2_used_leafs > m)
        {
            mid2.data += m;
            mid2.size = min(bi2_used_leafs - m, m);
            if (bi2_used_leafs > two_m)
            {
                high2.data += two_m;
                high2.size = bi2_used_leafs - two_m;
            }
            else
                high2.size = 0;
        }
        else
            mid2.size = high2.size = 0;

        /* for each quadratic polynomial (hx^2 + mx + l), using points (0, 1, -1, -2, infinity):
         *   p0 = l
         *   p1 = h + m + l
         *   p2 = h - m + l
         *   p3 = 4h - 2m + l
         *   p4 = h
         */

        /* calculate intermediate1 = h + l */
        if (high1.size)
            intermediate1 = bi_add(&high1, &low1);
        else
            intermediate1 = &low1;

        /* calculate intermediate2 = h + l */
        if (high2.size)
            intermediate2 = bi_add(&high2, &low2);
        else
            intermediate2 = &low2;

        /* calculate p0 = low1 * low2 */
        p0 = bi_toom_cook_mul(&low1, &low2, bi_used(&low1), bi_used(&low2));
        if (p0 == NULL) goto cleanup;

        /* calculate p1 = (low1 + mid1 + high1) * (low2 + mid2 + high2) */
        temp = bi_add(intermediate1, &mid1);
        if (temp == NULL) goto cleanup;

        temp2 = bi_add(intermediate2, &mid2);
        if (temp2 == NULL) goto cleanup;

        p1 = bi_toom_cook_mul(temp, temp2, bi_used(temp), bi_used(temp2));
        if (p1 == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;
        bi_destroy(temp2); temp2=NULL;

        /* calculate p2 = (low1 - mid1 + high1) * (low2 - mid2 + high2) */
        temp = bi_sub(intermediate1, &mid1);
        if (temp == NULL) goto cleanup;

        temp2 = bi_sub(intermediate2, &mid2);
        if (temp2 == NULL) goto cleanup;

        p2 = bi_toom_cook_mul(temp, temp2, bi_used(temp), bi_used(temp2));
        if (p2 == NULL) goto cleanup;

        if (intermediate1 != &low1) bi_destroy(intermediate1); intermediate1=NULL;
        if (intermediate2 != &low2) bi_destroy(intermediate2); intermediate2=NULL;
        // Do not cleanup temp and temp2 because the next step needs them

        /* calculate p3 = (low1 - 2*mid1 + 4*high1) * (low2 - 2*mid2 + 4*high2) */
        if ((temp = bi_add_assign(temp, &high1)) == NULL) goto cleanup;
        if ((temp = bi_shl1_assign(temp)) == NULL) goto cleanup;
        if ((temp = bi_sub_assign(temp, &low1)) == NULL) goto cleanup;

        if ((temp2 = bi_add_assign(temp2, &high2)) == NULL) goto cleanup;
        if ((temp2 = bi_shl1_assign(temp2)) == NULL) goto cleanup;
        if ((temp2 = bi_sub_assign(temp2, &low2)) == NULL) goto cleanup;

        p3 = bi_toom_cook_mul(temp, temp2, bi_used(temp), bi_used(temp2));
        if (p3 == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;
        bi_destroy(temp2); temp2=NULL;

        /* calculate p4 = high1 * high2 */
        p4 = bi_toom_cook_mul(&high1, &high2, bi_used(&high1), bi_used(&high2));
        if (p4 == NULL) goto cleanup;

        /* calculate interpolation */
        r0 = p0;
        r4 = p4;

        r3 = bi_sub(p3, p1);
        bi_destroy(p3); p3=NULL; // p3 no longer needed
        if (r3 == NULL) goto cleanup;
        if ((r3 = bi_div_3_assign(r3)) == NULL) goto cleanup;

        r1 = bi_sub(p1, p2);
        if (r1 == NULL) goto cleanup;
        (void) bi_shr_assign(r1, 1);
        bi_destroy(p1); p1=NULL; // p1 no longer needed

        r2 = bi_sub(p2, p0);
        if (r2 == NULL) goto cleanup;
        bi_destroy(p2); p2=NULL; // p2 no longer needed

        (void) bi_negate_assign(r3);
        if ((r3 = bi_add_assign(r3, r2)) == NULL) goto cleanup;
        (void) bi_shr_assign(r3, 1);

        temp = bi_shl1(r4);
        if ((r3 = bi_add_assign(r3, temp)) == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;

        if ((r2 = bi_add_assign(r2, r1)) == NULL) goto cleanup;
        if ((r2 = bi_sub_assign(r2, r4)) == NULL) goto cleanup;

        if ((r1 = bi_sub_assign(r1, r3)) == NULL) goto cleanup;

        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r3)) == NULL) goto cleanup;
        bi_destroy(r3); r3=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r2)) == NULL) goto cleanup;
        bi_destroy(r2); r2=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r1)) == NULL) goto cleanup;
        bi_destroy(r1); r1=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r0)) == NULL) goto cleanup;
        bi_destroy(r0); r0=NULL;

        // no need to free p0 because it's an alias for r0
        // no need to free p4 because we are returning it as r4 (r4 is an alias)
        // p1, p2, p3, temp, temp2, intermediate1, and intermediate2 are already freed
        r4->sign = bi->sign != bi2->sign;
        return r4;

cleanup:
        bi_destroy(r1);
        bi_destroy(r2);
        bi_destroy(r3);

        bi_destroy(p0);
        bi_destroy(temp);
        bi_destroy(temp2);
        if (intermediate1 != &low1) bi_destroy(intermediate1);
        if (intermediate2 != &low2) bi_destroy(intermediate2);
        bi_destroy(p1);
        bi_destroy(p2);
        bi_destroy(p3);
        bi_destroy(p4);
        return NULL;
    }
}
#endif

/* multiplies bi and bi2 and returns the result in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
bigint *bi_mul(const bigint *bi, const bigint *bi2)
{
    size_t bi_used_leafs = bi_used(bi), bi2_used_leafs = bi_used(bi2);
    bigint *result;

    /* handle multiply by 0, 1, or -1 */
    if (bi_used_leafs == 0 || bi2_used_leafs == 0) return bi_new();
    else if (bi_used_leafs == 1)
    {
        result = bi_mul_immediateu(bi2, bi->data[0]);
        if (result != NULL) result->sign = bi->sign != bi2->sign;
        return result;
    }
    else if (bi2_used_leafs == 1)
    {
        result = bi_mul_immediateu(bi, bi2->data[0]);
        if (result != NULL) result->sign = bi->sign != bi2->sign;
        return result;
    }

    if (bi_used_leafs >= BIGINT_TOOM_COOK_MIN_LEAFS && bi2_used_leafs >= BIGINT_TOOM_COOK_MIN_LEAFS)
#ifdef BIGINT_ENABLE_PTHREADS
        result = bi_toom_cook_mul(bi, bi2, bi_used_leafs, bi2_used_leafs, 1);
#else
        result = bi_toom_cook_mul(bi, bi2, bi_used_leafs, bi2_used_leafs);
#endif
    else if (bi_used_leafs >= BIGINT_KARATSUBA_MIN_LEAFS && bi2_used_leafs >= BIGINT_KARATSUBA_MIN_LEAFS)
#ifdef BIGINT_ENABLE_PTHREADS
        result = bi_karatsuba_mul(bi, bi2, bi_used_leafs, bi2_used_leafs, 1);
#else
        result = bi_karatsuba_mul(bi, bi2, bi_used_leafs, bi2_used_leafs);
#endif
    else
        result = bi_mul_internal(bi, bi2, bi_used_leafs, bi2_used_leafs);

    if (result == NULL) return NULL;
    result->sign = bi->sign != bi2->sign;
    return result;
}

/* TODO: can we multiply in place? */
/* multiplies bi and bi2 and assigns the product to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_mul_assign(bigint *bi, const bigint *bi2)
{
    bigint *result = bi_mul(bi, bi2);
    if (result != NULL)
    {
        bi_swap(bi, result);
        bi_destroy(result);
        return bi;
    }
    else
    {
        bi_destroy(bi);
        return NULL;
    }
}

/* multiplies bi and val and returns the result in bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_mul_immediateu_assign(bigint *bi, bi_leaf val)
{
    bi_dleaf carry = 0;
    size_t sz, i;

    if (val == 1) return bi;
    else if (val == 0) return bi_clear(bi);
    else if (ispow2x(val)) return bi_shl_assign(bi, lg2x(val));

    sz = bi_used(bi);

    /* handle multiply by 0, or multiply by 1/-1 */
    if (sz == 0) {bi_clear(bi); return bi;}
    else if (sz == 1 && bi->data[0] == 1) {bi->data[0] = val; return bi;}

    if (bi_resize(bi, sz + 1) == NULL) return NULL;

    for (i = 0; i < sz; ++i)
    {
        carry += (bi_dleaf) bi->data[i] * val;
        bi->data[i] = carry;
        carry >>= BIGINT_LEAFBITS;
    }

    bi->data[i] = carry;

    return bi;
}

/* multiplies bi and val and returns the product in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
bigint *bi_mul_immediateu(const bigint *bi, bi_leaf val)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_mul_immediateu_assign(result, val);
}

/* multiplies bi and val and returns the result in bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_mul_immediate_assign(bigint *bi, bi_signed_leaf val)
{
    bi_dleaf carry = 0;
    bi_leaf mul = val >= 0? (bi_leaf) val: -val == val? (bi_leaf) 1 << (sizeof(val)*CHAR_BIT - 1): (bi_leaf) -val;

    size_t sz, i;

    bi->sign = bi->sign != (val < 0);

    if (mul == 1) return bi;
    else if (mul == 0) return bi_clear(bi);

    sz = bi_used(bi);

    /* handle allocation failure, multiply by 0, or multiply by 1/-1 */
    if (sz == 0) return bi_clear(bi);
    else if (sz == 1 && bi->data[0] == 1) {bi->data[0] = mul; return bi;}

    if (bi_resize(bi, sz + 1) == NULL) return NULL;

    for (i = 0; i < sz; ++i)
    {
        carry += (bi_dleaf) bi->data[i] * mul;
        bi->data[i] = carry;
        carry >>= BIGINT_LEAFBITS;
    }

    bi->data[i] = carry;

    return bi;
}

/* multiplies bi and val and returns the product in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
bigint *bi_mul_immediate(const bigint *bi, bi_signed_leaf val)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_mul_immediate_assign(result, val);
}

static bigint *bi_square_internal(const bigint *bi, size_t used_leafs)
{
    bi_dleaf low_carry = 0, high_carry = 0;
    bigint *result;

    size_t i, j;
    size_t sz = used_leafs;

    if (sz == 0) return bi_new();
    else if (sz == 1)
    {
        result = bi_mul_immediateu(bi, bi->data[0]);
        if (result != NULL) result->sign = 0;
        return result;
    }

    result = bi_new_sized(max(sz * 2, BIGINT_MINLEAFS));
    if (result == NULL) return NULL;

    for (i = 0; i < sz * 2; ++i)
    {
        size_t j_start, j_end;

        j_start = i >= sz? i - sz + 1: 0;
        j_end = min(i, sz-1);

        for (j = j_start; j <= j_end; ++j)
        {
            bi_dleaf a, b, save;

            a = bi->data[j];
            b = bi->data[i - j];

            save = low_carry;
            low_carry += a * b;
            high_carry += low_carry < save;
        }

        result->data[i] = low_carry;
        low_carry = (low_carry >> BIGINT_LEAFBITS) | (high_carry << BIGINT_LEAFBITS);
        high_carry >>= BIGINT_LEAFBITS;
    }

    return result;
}

static bigint *bi_karatsuba_square(const bigint *bi, size_t bi_used_leafs);
static bigint *bi_toom_cook_square(const bigint *bi, size_t bi_used_leafs);

/* multiplies bi by itself and returns the result in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
bigint *bi_square(const bigint *bi)
{
    bigint *result;
    size_t bi_used_leafs = bi_used(bi);

    if (bi_used_leafs == 0) return bi_new();
    else if (bi_used_leafs == 1)
    {
        result = bi_mul_immediateu(bi, bi->data[0]);
        if (result != NULL) result->sign = 0;
        return result;
    }

    if (bi_used_leafs >= BIGINT_TOOM_COOK_MIN_LEAFS)
    {
        bigint b = *bi;
        b.sign = 0;

        result = bi_toom_cook_square(&b, bi_used_leafs);
    }
    else if (bi_used_leafs >= BIGINT_KARATSUBA_MIN_LEAFS)
    {
        bigint b = *bi;
        b.sign = 0;

        result = bi_karatsuba_square(&b, bi_used_leafs);
    }
    else
        result = bi_square_internal(bi, bi_used_leafs);

    if (result == NULL) return NULL;
    result->sign = 0;
    return result;
}

/* squares bi using Karatsuba multiplication, reverting to standard multiplication with small enough chunks */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_karatsuba_square(const bigint *bi, size_t bi_used_leafs)
{
    if (bi_used_leafs < BIGINT_KARATSUBA_MIN_LEAFS)
        return bi_square_internal(bi, bi_used_leafs);
    else
    {
        size_t m = bi_used_leafs / 2 + 1;
        bigint high, low, *z0 = NULL, *z1 = NULL, *z2 = NULL, *temp = NULL;

        low = high = *bi;
        low.size = m;
        high.data += m;
        high.size = bi_used_leafs - m;

        /* calculate z0 = low * low */
        z0 = bi_karatsuba_square(&low, bi_used(&low));
        if (z0 == NULL) goto cleanup;
        /* calculate z1 = (low + high) * (low + high) */
        temp = bi_add(&low, &high);
        if (temp == NULL) goto cleanup;
        z1 = bi_karatsuba_square(temp, bi_used(temp));
        if (z1 == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;
        /* calculate z2 = high * high */
        z2 = bi_karatsuba_square(&high, bi_used(&high));
        if (z2 == NULL) goto cleanup;
        /* calculate z1 = z1 - z2 - z0 */
        bi_fast_sub_assign(z1, z2);
        bi_fast_sub_assign(z1, z0);
        /* calculate z2 = z2 * 2^(2*m) + z1 * 2^m + z0 */
        if ((z2 = bi_shl_leaf_assign(z2, m)) == NULL) goto cleanup;
        if ((z2 = bi_add_assign(z2, z1)) == NULL) goto cleanup;
        bi_destroy(z1); z1=NULL;
        if ((z2 = bi_shl_leaf_assign(z2, m)) == NULL) goto cleanup;
        if ((z2 = bi_add_assign(z2, z0)) == NULL) goto cleanup;
        bi_destroy(z0);
        return z2;
cleanup:
        bi_destroy(z0);
        bi_destroy(temp);
        bi_destroy(z1);
        bi_destroy(z2);
        return NULL;
    }
}

/* squares bi using Toom-3 multiplication, reverting to standard multiplication with small enough chunks */
/* returns NULL on out of memory condition; the parameters are not destroyed */
static bigint *bi_toom_cook_square(const bigint *bi, size_t bi_used_leafs)
{
    if (bi_used_leafs < BIGINT_TOOM_COOK_MIN_LEAFS)
        return bi_karatsuba_square(bi, bi_used_leafs);
    else
    {
        size_t m = bi_used_leafs / 3 + 1;
        size_t two_m = 2 * m;
        bigint high, mid, low;
        bigint *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL, *p4 = NULL, *temp = NULL, *intermediate = NULL;
        bigint *r0 = NULL, *r1 = NULL, *r2 = NULL, *r3 = NULL, *r4 = NULL;

        low = mid = high = *bi;
        low.sign = mid.sign = high.sign = 0;
        low.size = min(m, bi_used_leafs);
        if (bi_used_leafs > m)
        {
            mid.data += m;
            mid.size = min(bi_used_leafs - m, m);
            if (bi_used_leafs > two_m)
            {
                high.data += two_m;
                high.size = bi_used_leafs - two_m;
            }
            else
                high.size = 0;
        }
        else
            mid.size = high.size = 0;

        /* for each quadratic polynomial (hx^2 + mx + l), using points (0, 1, -1, -2, infinity):
         *   p0 = l
         *   p1 = h + m + l
         *   p2 = h - m + l
         *   p3 = 4h - 2m + l
         *   p4 = h
         */

        /* calculate intermediate1 = h + l */
        if (high.size)
            intermediate = bi_add(&high, &low);
        else
            intermediate = &low;

        /* calculate p0 = low1 * low2 */
        p0 = bi_toom_cook_square(&low, bi_used(&low));
        if (p0 == NULL) goto cleanup;

        /* calculate p1 = (low1 + mid1 + high1) * (low2 + mid2 + high2) */
        temp = bi_add(intermediate, &mid);
        if (temp == NULL) goto cleanup;

        p1 = bi_toom_cook_square(temp, bi_used(temp));
        if (p1 == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;

        /* calculate p2 = (low1 - mid1 + high1) * (low2 - mid2 + high2) */
        temp = bi_sub(intermediate, &mid);
        if (temp == NULL) goto cleanup;

        p2 = bi_toom_cook_square(temp, bi_used(temp));
        if (p2 == NULL) goto cleanup;

        if (intermediate != &low) bi_destroy(intermediate); intermediate=NULL;
        // Do not cleanup temp because the next step needs it

        /* calculate p3 = (low1 - 2*mid1 + 4*high1) * (low2 - 2*mid2 + 4*high2) */
        if ((temp = bi_add_assign(temp, &high)) == NULL) goto cleanup;
        if ((temp = bi_shl1_assign(temp)) == NULL) goto cleanup;
        if ((temp = bi_sub_assign(temp, &low)) == NULL) goto cleanup;

        p3 = bi_toom_cook_square(temp, bi_used(temp));
        if (p3 == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;

        /* calculate p4 = high1 * high2 */
        p4 = bi_toom_cook_square(&high, bi_used(&high));
        if (p4 == NULL) goto cleanup;

        /* calculate interpolation */
        r0 = p0;
        r4 = p4;

        r3 = bi_sub(p3, p1);
        bi_destroy(p3); p3=NULL; // p3 no longer needed
        if (r3 == NULL) goto cleanup;
        // TODO: speed up division by 3?
        if ((r3 = bi_div_3_assign(r3)) == NULL) goto cleanup;

        r1 = bi_sub(p1, p2);
        if (r1 == NULL) goto cleanup;
        (void) bi_shr_assign(r1, 1);
        bi_destroy(p1); p1=NULL; // p1 no longer needed

        r2 = bi_sub(p2, p0);
        if (r2 == NULL) goto cleanup;
        bi_destroy(p2); p2=NULL; // p2 no longer needed

        (void) bi_negate_assign(r3);
        if ((r3 = bi_add_assign(r3, r2)) == NULL) goto cleanup;
        (void) bi_shr_assign(r3, 1);

        temp = bi_shl1(r4);
        if ((r3 = bi_add_assign(r3, temp)) == NULL) goto cleanup;
        bi_destroy(temp); temp=NULL;

        if ((r2 = bi_add_assign(r2, r1)) == NULL) goto cleanup;
        if ((r2 = bi_sub_assign(r2, r4)) == NULL) goto cleanup;

        if ((r1 = bi_sub_assign(r1, r3)) == NULL) goto cleanup;

        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r3)) == NULL) goto cleanup;
        bi_destroy(r3); r3=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r2)) == NULL) goto cleanup;
        bi_destroy(r2); r2=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r1)) == NULL) goto cleanup;
        bi_destroy(r1); r1=NULL;
        if ((r4 = bi_shl_leaf_assign(r4, m)) == NULL) goto cleanup;
        if ((r4 = bi_add_assign(r4, r0)) == NULL) goto cleanup;
        bi_destroy(r0); r0=NULL;

        // no need to free p0 because it's an alias for r0
        // no need to free p4 because we are returning it as r4 (r4 is an alias)
        // p1, p2, p3, temp, temp2, intermediate1, and intermediate2 are already freed
        r4->sign = 0;
        return r4;

cleanup:
        bi_destroy(r1);
        bi_destroy(r2);
        bi_destroy(r3);

        bi_destroy(p0);
        bi_destroy(temp);
        if (intermediate != &low) bi_destroy(intermediate);
        bi_destroy(p1);
        bi_destroy(p2);
        bi_destroy(p3);
        bi_destroy(p4);
        return NULL;
    }
}

/* TODO: can we multiply in place? */
/* multiplies bi by itself and assigns the product to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_square_assign(bigint *bi)
{
    bigint *result = bi_square(bi);
    if (result != NULL)
    {
        bi_swap(bi, result);
        bi_destroy(result);
        return bi;
    }
    else
    {
        bi_destroy(bi);
        return NULL;
    }
}

/* calculates square root of bi and returns the result in a new bigint */
/* returns NULL on out of memory condition; the parameters are not destroyed */
bigint *bi_sqrt(const bigint *bi)
{
    size_t sz = bi_bitcount(bi);
    bigint *val = bi_copy(bi);
    bigint *root = bi_new();
    bigint *bit;

    if (sz == 0 || bi_is_negative(bi))
        return root;

    bit = bi_new();
    if (root == NULL || bit == NULL)
        goto cleanup;

    bi_assign(bit, 1);
    if ((bit = bi_shl_assign(bit, (sz-1)/2*2)) == NULL) goto cleanup;

    while (bi_cmp(bit, bi) > 0)
        if ((bit = bi_shr_assign(bit, 2)) == NULL) goto cleanup;

    while (!bi_is_zero(bit))
    {
        bigint *temp = bi_add(root, bit);
        if (temp == NULL) goto cleanup;

        if (bi_cmp_mag(val, temp) >= 0)
        {
            if ((val = bi_sub_assign(val, temp)) == NULL ||
                (root = bi_shr_assign(root, 1)) == NULL ||
                (root = bi_add_assign(root, bit)) == NULL) goto cleanup;
        }
        else
            if ((root = bi_shr_assign(root, 1)) == NULL) goto cleanup;
        bi_destroy(temp);
        if ((bit = bi_shr_assign(bit, 2)) == NULL) goto cleanup;
    }

    bi_destroy(val);
    bi_destroy(bit);
    return root;

cleanup:
    bi_destroy(root);
    bi_destroy(val);
    bi_destroy(bit);
    return NULL;
}

/* calculates square root of bi and assigns the product to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_sqrt_assign(bigint *bi)
{
    bigint *result = bi_sqrt(bi);
    if (result != NULL)
    {
        bi_swap(bi, result);
        bi_destroy(result);
        return bi;
    }
    else
    {
        bi_destroy(bi);
        return NULL;
    }
}

/* calculates factorial of n and returns the result in a new bigint */
/* returns NULL on out of memory condition */
bigint *bi_fact(bi_uintmax n)
{
    bigint *result = bi_new(), *mult = bi_new();
    bi_uintmax i = 1;

    if (result == NULL || mult == NULL) {bi_destroy(result); bi_destroy(mult); return NULL;}
    bi_assign(result, 1);

    for (; i < n; ++i)
    {
        if (bi_assignlu(mult, i+1) == NULL) {bi_destroy(result); return NULL;}
        if (bi_mul_assign(result, mult) == NULL) {bi_destroy(mult); return NULL;}
    }

    return result;
}

/* raises bi to the nth power and returns the result in a new bigint */
/* returns NULL on out of memory condition */
bigint *bi_uexp(const bigint *bi, bi_uintmax n)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_uexp_assign(result, n);
}

/* raises bi to the nth power and assigns the result to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_uexp_assign(bigint *bi, bi_uintmax n)
{
    bigint *y;

    if (n == 1) return bi;
    else if (n == 2) return bi_square_assign(bi);
    else if (n == 0) return bi_assignu(bi, 1);
    else if (bi_is_one(bi)) return bi;

    y = bi_new();
    if (y == NULL || bi_assign(y, 1) == NULL)
    {
        bi_destroy(y);
        return NULL;
    }

    while (n > 1)
    {
        if ((n & 1) && bi_mul_assign(y, bi) == NULL) {bi_destroy(bi); return NULL;}
        n >>= 1;
        if (bi_square_assign(bi) == NULL) {bi_destroy(y); return NULL;}
    }

    if (bi_mul_assign(bi, y) == NULL) {bi_destroy(y); return NULL;}
    bi_destroy(y);
    return bi;
}

/* raises bi to the nth power and returns the result in a new bigint */
/* returns NULL on out of memory condition */
bigint *bi_exp(const bigint *bi, bi_intmax n)
{
    if (n < 0) return bi_new();
    return bi_uexp(bi, n);
}

/* raises bi to the nth power and assigns the result to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_exp_assign(bigint *bi, bi_intmax n)
{
    if (n < 0) return bi_clear(bi);
    return bi_uexp_assign(bi, n);
}

/* raises bi to the nth power and returns the result (using modulo `mod`) in a new bigint */
/* returns NULL on out of memory condition */
bigint *bi_uexp_mod(const bigint *bi, bi_uintmax n, const bigint *mod)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_uexp_mod_assign(result, n, mod);
}

/* raises bi to the nth power and assigns the result (using modulo `mod`) to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_uexp_mod_assign(bigint *bi, bi_uintmax n, const bigint *mod)
{
    bigint *y;

    if (bi_is_one(mod)) return bi_clear(bi);

    if (n == 1) return bi_mod_assign(bi, mod);
    else if (n == 2) return bi_square_assign(bi) == NULL? NULL: bi_mod_assign(bi, mod);
    else if (n == 0) return bi_assignu(bi, 1);
    else if (bi_is_one(bi)) return bi;

    y = bi_new();
    if (y == NULL || bi_assign(y, 1) == NULL)
    {
        bi_destroy(y);
        return NULL;
    }

    if (bi_mod_assign(bi, mod) == NULL) {bi_destroy(y); return NULL;}

    while (n > 0)
    {
        if ((n & 1) &&
            (bi_mul_assign(y, bi) == NULL ||
             bi_mod_assign(y, mod) == NULL)) {bi_destroy(bi); return NULL;}
        n >>= 1;
        if (bi_square_assign(bi) == NULL ||
            bi_mod_assign(bi, mod) == NULL) {bi_destroy(y); return NULL;}
    }

    bi_swap(bi, y);
    bi_destroy(y);
    return bi;
}

/* raises bi to the nth power and returns the result (using modulo `mod`) in a new bigint */
/* returns NULL on out of memory condition */
bigint *bi_exp_mod(const bigint *bi, bi_intmax n, const bigint *mod)
{
    bigint *result;
    if (n < 0) return bi_new();
    if ((result = bi_copy(bi)) == NULL) return NULL;
    return bi_exp_mod_assign(result, n, mod);
}

/* raises bi to the nth power and assigns the result (using modulo `mod`) to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_exp_mod_assign(bigint *bi, bi_intmax n, const bigint *mod)
{
    bigint *result;
    if (n < 0) return bi_clear(bi);
    if ((result = bi_copy(bi)) == NULL) return NULL;
    return bi_uexp_mod_assign(result, n, mod);
}

/* raises bi to the nth power and returns the result (using modulo `mod`) in a new bigint */
/* returns NULL on out of memory condition */
bigint *bi_large_exp_mod(const bigint *bi, const bigint *n, const bigint *mod)
{
    bigint *result = bi_copy(bi);
    if (result == NULL) return NULL;
    return bi_large_exp_mod_assign(result, n, mod);
}

// TODO: not yet implemented
/* raises bi to the nth power and assigns the result (using modulo `mod`) to bi, returns bi */
/* returns NULL and destroys bi on out of memory condition */
bigint *bi_large_exp_mod_assign(bigint *bi, const bigint *n, const bigint *mod)
{
    bigint *y;
    size_t bits, i;

    if (bi_is_negative(n) || bi_is_one(mod)) return bi_clear(bi);

    if (bi_is_one(n)) return bi_mod_assign(bi, mod);
    else if (bi_cmp_immu(n, 2) == 0) return bi_square_assign(bi) == NULL? NULL: bi_mod_assign(bi, mod);
    else if (bi_is_zero(n)) return bi_assignu(bi, 1);
    else if (bi_is_one(bi)) return bi;

    y = bi_new();
    if (y == NULL || bi_assign(y, 1) == NULL)
    {
        bi_destroy(y);
        return NULL;
    }

    if (bi_mod_assign(bi, mod) == NULL) {bi_destroy(y); return NULL;}

    bits = bi_bitcount(n);
    for (i = 0; i < bits; ++i)
    {
        if (bi_bit(n, i) &&
            (bi_mul_assign(y, bi) == NULL ||
             bi_mod_assign(y, mod) == NULL)) {bi_destroy(bi); return NULL;}

        if (bi_square_assign(bi) == NULL ||
            bi_mod_assign(bi, mod) == NULL) {bi_destroy(y); return NULL;}
    }

    bi_swap(bi, y);
    bi_destroy(y);
    return bi;
}

static bigint *bi_divmod_immediate_internal(bigint *dst, const bigint *bi, bi_signed_leaf denom, bigint **q, bi_signed_leaf *r);

/* divides bi by bi2 and sets *q to a new bigint containing the quotient, *r to a new bigint containing the remainder, returns quotient */
/* if dst == NULL, the quotient is put in a new bigint, otherwise, the quotient is placed in dst (without testing size of dst) */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed and q and r are not set */
static bigint *bi_divmod_internal(bigint *dst, const bigint *bi, const bigint *bi2, bigint **q, bigint **r)
{
    bigint *quotient;
    bigint *remainder;
    size_t i, quotient_zero = 1;

    i = bi_used(bi2);
    if (i <= 1)
    {
        if (i == 0) return NULL;
        if ((bi2->data[0] & BIGINT_SIGNED_LEAFMAX) == bi2->data[0])
        {
            bi_signed_leaf l;
            if (bi_divmod_immediate_internal(dst, bi, bi_to_int(bi2), q, &l) == NULL) return NULL;
            *r = bi_new();
            if (*r == NULL) {bi_destroy(*q); return NULL;}
            bi_assign(*r, l);
            return *q;
        }
    }

    i = bi_bitcount(bi);
    if (dst == NULL)
    {
        quotient = bi_new_sized(max((i+BIGINT_LEAFBITS-1) / BIGINT_LEAFBITS, BIGINT_MINLEAFS));
        if (quotient == NULL) return NULL;
    }
    else
        quotient = dst;

    remainder = bi_new();
    if (remainder == NULL) goto cleanup;

    for (; i > 0; --i)
    {
        size_t byte_index = (i-1) / BIGINT_LEAFBITS;
        size_t bit_index = (i-1) % BIGINT_LEAFBITS;

        if ((remainder = bi_shl1_assign(remainder)) == NULL) goto cleanup;
        remainder->data[0] |= (bi->data[byte_index] >> bit_index) & 1;
        if (bi_cmp_mag(remainder, bi2) >= 0)
        {
            bi_fast_sub_assign(remainder, bi2);
            quotient->data[byte_index] |= (bi_leaf) 1 << bit_index;
            quotient_zero = 0;
        }
        else
            quotient->data[byte_index] &= ~((bi_leaf) 1 << bit_index);
    }

    remainder->sign = bi->sign && !bi_is_zero(remainder);
    quotient->sign = bi->sign != bi2->sign && !quotient_zero;
    *q = quotient;
    *r = remainder;
    return quotient;

cleanup:
    bi_destroy(remainder);
    bi_destroy(quotient);
    return NULL;
}

/* divides bi by bi2 and sets *q to a new bigint containing the quotient, *r to a new bigint containing the remainder, returns quotient */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed and q and r are not set */
bigint *bi_divmod(const bigint *bi, const bigint *bi2, bigint **q, bigint **r)
{
    return bi_divmod_internal(NULL, bi, bi2, q, r);
}

/* divides bi by bi2 and places the quotient in bi and the remainder in r, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_divmod_assign(bigint *bi, const bigint *bi2, bigint **r)
{
    bigint *q;
    return bi_divmod_internal(bi, bi, bi2, &q, r);
}

/* divides bi by bi2 and returns the quotient in a new bigint */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed */
bigint *bi_div(const bigint *bi, const bigint *bi2)
{
    bigint *q, *r;
    if (bi_divmod(bi, bi2, &q, &r) == NULL) return NULL;
    bi_destroy(r);
    return q;
}

/* divides bi by bi2 and places the quotient in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_div_assign(bigint *bi, const bigint *bi2)
{
    bigint *q, *r;
    if (bi_divmod_internal(bi, bi, bi2, &q, &r) == NULL) return NULL;
    bi_destroy(r);
    return bi;
}

/* divides bi by bi2 and returns the remainder in a new bigint */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed */
bigint *bi_mod(const bigint *bi, const bigint *bi2)
{
    bigint *q, *r;
    if (bi_divmod(bi, bi2, &q, &r) == NULL) return NULL;
    bi_destroy(q);
    return r;
}

/* TODO: can we modulo in place? */
/* divides bi by bi2 and places the remainder in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_mod_assign(bigint *bi, const bigint *bi2)
{
    bigint *q, *r;
    if (bi_divmod(bi, bi2, &q, &r) == NULL)
    {
        bi_destroy(bi);
        return NULL;
    }
    bi_destroy(q);
    bi_swap(bi, r);
    bi_destroy(r);
    return bi;
}

/* divides bi by 3 and sets *q to a new bigint containing the quotient, *r to the remainder, returns quotient */
/* if dst == NULL, the quotient is put in a new bigint, otherwise, the quotient is placed in dst (without testing size of dst) */
/* returns NULL on out of memory; the parameters are not destroyed and q and r are not set */
static bigint *bi_divmod_3_internal(bigint *dst, const bigint *bi, bigint **q, bi_signed_leaf *r)
{
    bigint *quotient;
    size_t i;
    bi_leaf remainder = 0;

    i = bi_used(bi);
    if (dst == NULL)
    {
        quotient = bi_new_sized(max(i, BIGINT_MINLEAFS));
        if (quotient == NULL) return NULL;
    }
    else
        quotient = dst;

    for (; i > 0; --i)
    {
        bi_dleaf numer = remainder;
        numer <<= BIGINT_LEAFBITS;
        numer |= bi->data[i-1];

        // Code taken from http://www.hackersdelight.org/hdcodetxt/divuc.c.txt
#if CHAR_BIT == 8
        {
            bi_dleaf q = numer;
            q = (numer >> 2) + (numer >> 4);
            q += q >> 4;
#if BIGINT_LEAFBITS*2 > 8
            q += q >> 8;
#if BIGINT_LEAFBITS*2 > 16
            q += q >> 16;
#if BIGINT_LEAFBITS*2 > 32
            q += q >> 32;
#if BIGINT_LEAFBITS*2 > 64
            q += q >> 64;
#endif
#endif
#endif
#endif
            remainder = numer - (q + (q << 1));
            quotient->data[i-1] = q = q + ((remainder + 5 + (remainder << 2)) >> 4);
            remainder = numer - (q + (q << 1));
        }
#else
        quotient->data[i-1] = denom / 3;
        remainder = denom % 3;
#endif
    }

    if (bi->sign)
        *r = -((bi_signed_leaf) remainder);
    else
        *r = remainder;
    quotient->sign = bi->sign && !bi_is_zero(quotient);
    *q = quotient;
    return quotient;
}

/* divides bi by 3 and sets *q to a new bigint containing the quotient, *r to the remainder, returns quotient */
/* returns NULL on out of memory; the parameters are not destroyed and q and r are not set */
bigint *bi_divmod_3(const bigint *bi, bigint **q, bi_signed_leaf *r)
{
    return bi_divmod_3_internal(NULL, bi, q, r);
}

/* divides bi by 3 and returns the quotient in a new bigint */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed */
bigint *bi_div_3(const bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    return bi_divmod_3(bi, &q, &r);
}

/* divides bi by 3 and places the quotient in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_div_3_assign(bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    return bi_divmod_3_internal(bi, bi, &q, &r);
}

/* divides bi by 3 and places the quotient in bi and the remainder in r, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_divmod_3_assign(bigint *bi, bi_signed_leaf *r)
{
    bigint *q;
    return bi_divmod_3_internal(bi, bi, &q, r);
}

/* divides bi by 3 and returns the remainder in a new bigint */
/* returns 0 on out of memory or divide by zero; the parameters are not destroyed */
bi_signed_leaf bi_mod_3(const bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    if (bi_divmod_3(bi, &q, &r) == NULL) return 0;
    bi_destroy(q);
    return r;
}

/* TODO: can we modulo in place? */
/* divides bi by 3 and places the remainder in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_mod_3_assign(bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    if (bi_divmod_3(bi, &q, &r) == NULL)
    {
        bi_destroy(bi);
        return NULL;
    }
    bi_destroy(q);
    return bi_assign(bi, r);
}

/* divides bi by 10 and sets *q to a new bigint containing the quotient, *r to the remainder, returns quotient */
/* if dst == NULL, the quotient is put in a new bigint, otherwise, the quotient is placed in dst (without testing size of dst) */
/* returns NULL on out of memory; the parameters are not destroyed and q and r are not set */
static bigint *bi_divmod_10_internal(bigint *dst, const bigint *bi, bigint **q, bi_signed_leaf *r)
{
    bigint *quotient;
    size_t i;
    bi_leaf remainder = 0;

    i = bi_used(bi);
    if (dst == NULL)
    {
        quotient = bi_new_sized(max(i, BIGINT_MINLEAFS));
        if (quotient == NULL) return NULL;
    }
    else
        quotient = dst;

    for (; i > 0; --i)
    {
        bi_dleaf numer = remainder;
        numer <<= BIGINT_LEAFBITS;
        numer |= bi->data[i-1];

        // Code taken from http://www.hackersdelight.org/hdcodetxt/divuc.c.txt
#if CHAR_BIT == 8
        {
            bi_dleaf q = numer;
            q = (numer >> 1) + (numer >> 2);
            q += q >> 4;
#if BIGINT_LEAFBITS*2 > 8
            q += q >> 8;
#if BIGINT_LEAFBITS*2 > 16
            q += q >> 16;
#if BIGINT_LEAFBITS*2 > 32
            q += q >> 32;
#if BIGINT_LEAFBITS*2 > 64
            q += q >> 64;
#endif
#endif
#endif
#endif
            q >>= 3;
            remainder = numer - ((q + (q << 2)) << 1);
            quotient->data[i-1] = q = q + ((remainder + 6) >> 4);
            remainder = numer - ((q + (q << 2)) << 1);
        }
#else
        quotient->data[i-1] = numer / 10;
        remainder = numer % 10;
#endif
    }

    if (bi->sign)
        *r = -((bi_signed_leaf) remainder);
    else
        *r = remainder;
    quotient->sign = bi->sign && !bi_is_zero(quotient);
    *q = quotient;
    return quotient;
}

/* divides bi by 10 and sets *q to a new bigint containing the quotient, *r to the remainder, returns quotient */
/* returns NULL on out of memory; the parameters are not destroyed and q and r are not set */
bigint *bi_divmod_10(const bigint *bi, bigint **q, bi_signed_leaf *r)
{
    return bi_divmod_10_internal(NULL, bi, q, r);
}

/* divides bi by 10 and returns the quotient in a new bigint */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed */
bigint *bi_div_10(const bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    return bi_divmod_10(bi, &q, &r);
}

/* divides bi by 10 and places the quotient in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_div_10_assign(bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    return bi_divmod_10_internal(bi, bi, &q, &r);
}

/* divides bi by 10 and places the quotient in bi and the remainder in r, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_divmod_10_assign(bigint *bi, bi_signed_leaf *r)
{
    bigint *q;
    return bi_divmod_10_internal(bi, bi, &q, r);
}

/* divides bi by 10 and returns the remainder in a new bigint */
/* returns 0 on out of memory or divide by zero; the parameters are not destroyed */
bi_signed_leaf bi_mod_10(const bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    if (bi_divmod_10(bi, &q, &r) == NULL) return 0;
    bi_destroy(q);
    return r;
}

/* TODO: can we modulo in place? */
/* divides bi by 10 and places the remainder in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_mod_10_assign(bigint *bi)
{
    bigint *q;
    bi_signed_leaf r;
    if (bi_divmod_10(bi, &q, &r) == NULL)
    {
        bi_destroy(bi);
        return NULL;
    }
    bi_destroy(q);
    return bi_assign(bi, r);
}

/* divides bi by denom and sets *q to a new bigint containing the quotient, *r to the remainder, returns quotient */
/* if dst == NULL, the quotient is put in a new bigint, otherwise, the quotient is placed in dst (without testing size of dst) */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed and q and r are not set */
static bigint *bi_divmod_immediate_internal(bigint *dst, const bigint *bi, bi_signed_leaf denom, bigint **q, bi_signed_leaf *r)
{
    bigint *quotient;
    size_t i;
    bi_leaf remainder = 0, d = denom >= 0? (bi_leaf) denom: -denom == denom? (bi_leaf) 1 << (sizeof(denom)*CHAR_BIT - 1): (bi_leaf) -denom;

    if (d == 0) return NULL;
    else if (d == 1)
    {
        *r = 0;
        if (dst == NULL)
            return *q = denom == 1? bi_copy(bi): bi_negate(bi);
        else if (bi_copy_to(dst, bi) != NULL)
            return *q = denom == 1? dst: bi_negate_assign(dst);
        else
            return NULL;
    }
    else if (ispow2x(d))
    {
        *r = bi->sign? -(bi_signed_leaf) (bi_to_intu(bi) & (d-1)): (bi_signed_leaf) (bi_to_intu(bi) & (d-1));
        if (dst == NULL)
            *q = bi_shr(bi, lg2x(d));
        else if (bi_copy_to(dst, bi) != NULL)
            *q = bi_shr_assign(dst, lg2x(d));
        else
            return *q = NULL;

        if (*q != NULL)
            (*q)->sign = bi->sign != (denom < 0) && !bi_is_zero(*q);
        return *q;
    }
    else if (denom == 3)
        return bi_divmod_3_internal(dst, bi, q, r);
    else if (denom == 10)
        return bi_divmod_10_internal(dst, bi, q, r);

#if 0
    /* bit by bit version */
    i = bi_bitcount(bi);
    /* shortcut for single-leaf division */
    if (i <= BIGINT_LEAFBITS)
    {
        bi_leaf n = bi->data[0];
        if (dst == NULL)
        {
            quotient = bi_new();
            if (quotient == NULL) return NULL;
        }
        else
            quotient = dst;
        quotient->data[0] = n / d;
        remainder = n % d;
    }
    else
    {
        if (dst == NULL)
        {
            quotient = bi_new_sized(max((i+BIGINT_LEAFBITS-1) / BIGINT_LEAFBITS, BIGINT_MINLEAFS));
            if (quotient == NULL) return NULL;
        }
        else
            quotient = dst;

        for (; i > 0; --i)
        {
            size_t byte_index = (i-1) / BIGINT_LEAFBITS;
            size_t bit_index = (i-1) % BIGINT_LEAFBITS;

            remainder <<= 1;
            remainder |= (bi->data[byte_index] >> bit_index) & 1;
            if (remainder >= d)
            {
                remainder -= d;
                quotient->data[byte_index] |= (bi_leaf) 1 << bit_index;
            }
            else
                quotient->data[byte_index] &= ~((bi_leaf) 1 << bit_index);
        }
    }
#else
    /* digit by digit (leaf by leaf) version */
    i = bi_used(bi);
    if (dst == NULL)
    {
        quotient = bi_new_sized(max(i, BIGINT_MINLEAFS));
        if (quotient == NULL) return NULL;
    }
    else
        quotient = dst;

    for (; i > 0; --i)
    {
        bi_dleaf numer = remainder;
        numer <<= BIGINT_LEAFBITS;
        numer |= bi->data[i-1];

        quotient->data[i-1] = numer / d;
        remainder = numer % d;
    }
#endif
    if (bi->sign)
        *r = -((bi_signed_leaf) remainder);
    else
        *r = remainder;
    quotient->sign = bi->sign != (denom < 0) && !bi_is_zero(quotient);
    *q = quotient;
    return quotient;
}

/* divides bi by denom and sets *q to a new bigint containing the quotient, *r to the remainder, returns quotient */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed and q and r are not set */
bigint *bi_divmod_immediate(const bigint *bi, bi_signed_leaf denom, bigint **q, bi_signed_leaf *r)
{
    return bi_divmod_immediate_internal(NULL, bi, denom, q, r);
}

/* divides bi by denom and returns the quotient in a new bigint */
/* returns NULL on out of memory or divide by zero; the parameters are not destroyed */
bigint *bi_div_immediate(const bigint *bi, bi_signed_leaf denom)
{
    bigint *q;
    bi_signed_leaf r;
    return bi_divmod_immediate(bi, denom, &q, &r);
}

/* divides bi by denom and places the quotient in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_div_immediate_assign(bigint *bi, bi_signed_leaf denom)
{
    bigint *q;
    bi_signed_leaf r;
    return bi_divmod_immediate_internal(bi, bi, denom, &q, &r);
}

/* divides bi by denom and places the quotient in bi and the remainder in r, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_divmod_immediate_assign(bigint *bi, bi_signed_leaf denom, bi_signed_leaf *r)
{
    bigint *q;
    return bi_divmod_immediate_internal(bi, bi, denom, &q, r);
}

/* divides bi by denom and returns the remainder in a new bigint */
/* returns 0 on out of memory or divide by zero; the parameters are not destroyed */
bi_signed_leaf bi_mod_immediate(const bigint *bi, bi_signed_leaf denom)
{
    bigint *q;
    bi_signed_leaf r;
    if (bi_divmod_immediate(bi, denom, &q, &r) == NULL) return 0;
    bi_destroy(q);
    return r;
}

/* TODO: can we modulo in place? */
/* divides bi by denom and places the remainder in bi, returns bi */
/* returns NULL and destroys bi on out of memory or divide by zero */
bigint *bi_mod_immediate_assign(bigint *bi, bi_signed_leaf denom)
{
    bigint *q;
    bi_signed_leaf r;
    if (bi_divmod_immediate(bi, denom, &q, &r) == NULL)
    {
        bi_destroy(bi);
        return NULL;
    }
    bi_destroy(q);
    return bi_assign(bi, r);
}

/* finds greatest common divisor of absolute values of bi and bi2 */
/* returns NULL on out of memory */
/* algorithm taken from https://en.wikipedia.org/wiki/Binary_GCD_algorithm */
bigint *bi_gcd(const bigint *bi, const bigint *bi2)
{
    bigint *u = NULL, *v = NULL;
    size_t pow = 0;

    if (bi_cmp(bi, bi2) == 0 || bi_is_zero(bi2))
        return bi_copy(bi);
    else if (bi_is_zero(bi))
        return bi_copy(bi2);

    u = bi_abs(bi);
    v = bi_abs(bi2);
    if (u == NULL || v == NULL) goto cleanup;

    // extract even factor to pow
    while ((u->data[0] & 1) == 0 && (v->data[0] & 1) == 0)
    {
        ++pow;
        u = bi_shr_assign(u, 1);
        v = bi_shr_assign(v, 1);
        if (u == NULL || v == NULL) goto cleanup;
    }

    // remove even factors of u
    while ((u->data[0] & 1) == 0)
    {
        u = bi_shr_assign(u, 1);
        if (u == NULL) goto cleanup;
    }

    do
    {
        while ((v->data[0] & 1) == 0)
        {
            v = bi_shr_assign(v, 1);
            if (v == NULL) goto cleanup;
        }

        if (bi_cmp_mag(u, v) > 0)
            bi_swap(u, v);

        v = bi_fast_sub_assign(v, u);
        if (v == NULL) goto cleanup;
    } while (!bi_is_zero(v));

    u = bi_shl_assign(u, pow);
    if (u != NULL) return u;

cleanup:
    bi_destroy(u);
    bi_destroy(v);
    return NULL;
}

/* swaps contents of bigints a and b */
void bi_swap(bigint *bi_a, bigint *bi_b)
{
    bigint temp = *bi_a;
    *bi_a = *bi_b;
    *bi_b = temp;
}

static void bi_fprint_basepow2(FILE *f, const bigint *bi, char *array, size_t array_size, size_t lg2_base)
{
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    bigint *cpy = NULL;
    bi_leaf mask = ((bi_leaf) 1 << lg2_base) - 1;
    size_t size = array_size, ptr = 0;

    if (mask == 7 || (mask >= 31 && mask != 255))
    {
        cpy = bi_copy(bi);
        if (cpy == NULL) goto cleanup;
        while (!bi_is_zero(cpy))
        {
            if (ptr == size)
            {
                void *temp;
                size += size >> 1;
                temp = realloc(array, size);
                if (temp == NULL)
                    goto cleanup;
                array = temp;
            }
            array[ptr++] = mask > 31? (char) (cpy->data[0] & mask): alphabet[cpy->data[0] & mask];
            if ((cpy = bi_shr_assign(cpy, lg2_base)) == NULL) goto cleanup;
        }

        while (ptr > 0)
            fputc(array[--ptr], f);
    }
    else
    {
        size_t i = bi->size, j;
        int started = 0;
        for (; i > 0; --i)
        {
            for (j = BIGINT_LEAFBITS; j > 0; j -= lg2_base)
            {
                int data = (bi->data[i-1] >> (j-lg2_base)) & mask;
                started |= data;
                if (started)
                    fputc(mask == 255? data: alphabet[data], f);
            }
        }
    }

cleanup:
    free(array);
    bi_destroy(cpy);
}

/* prints bi to stream f, with specified base */
void bi_fprint(FILE *f, const bigint *bi, size_t base)
{
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    bigint *cpy = NULL;
    size_t size = 100, ptr = 0;
    bi_signed_leaf r;
    char *array;

    array = malloc(size);
    if (array == NULL) goto cleanup;

    if (bi->sign)
        fputc('-', f);

    if (bi_is_zero(bi))
        fputc('0', f);
    else if (base < 2)
        ; // do nothing
    else if (ispow2(base))
    {
        bi_fprint_basepow2(f, bi, array, size, lg2(base));
        cpy = NULL;
        array = NULL;
    }
    else
    {
        cpy = bi_copy(bi);
        if (cpy == NULL) goto cleanup;
        cpy->sign = 0;

        while (!bi_is_zero(cpy))
        {
            if (ptr == size)
            {
                void *temp;
                size += size >> 1;
                temp = realloc(array, size);
                if (temp == NULL)
                    goto cleanup;
                array = temp;
            }
            if ((cpy = bi_divmod_immediate_assign(cpy, base, &r)) == NULL) goto cleanup;
            array[ptr++] = r <= 36? alphabet[r]: r;
        }

        while (ptr > 0)
            fputc(array[--ptr], f);
    }

cleanup:
    free(array);
    bi_destroy(cpy);
}

static bigint_string bi_sprint_basepow2(const bigint *bi, char *array, size_t array_start, size_t array_size, size_t lg2_base)
{
    char *result = NULL;
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    bigint *cpy = NULL;
    bi_leaf mask = ((bi_leaf) 1 << lg2_base) - 1;
    size_t size = array_size, ptr = array_start;

    if (mask == 7 || (mask >= 31 && mask != 255))
    {
        cpy = bi_copy(bi);
        if (cpy == NULL) goto cleanup;
        while (!bi_is_zero(cpy))
        {
            if (ptr == size)
            {
                void *temp;
                size += size >> 1;
                temp = realloc(array, size);
                if (temp == NULL)
                    goto cleanup;
                array = temp;
            }
            array[ptr++] = mask > 31? (char) (cpy->data[0] & mask): alphabet[cpy->data[0] & mask];
            if ((cpy = bi_shr_assign(cpy, lg2_base)) == NULL) goto cleanup;
        }

        result = malloc(size+1);
        if (result == NULL)
            goto cleanup;

        size = 0;
        while (ptr > 0)
            result[size++] = array[--ptr];
        result[size] = 0;
        ptr = size;
    }
    else
    {
        size_t i = bi->size, j;
        int started = 0;
        for (; i > 0; --i)
        {
            for (j = BIGINT_LEAFBITS; j > 0; j -= lg2_base)
            {
                int data = (bi->data[i-1] >> (j-lg2_base)) & mask;
                started |= data;
                if (started)
                {
                    if (ptr + 1 >= size)
                    {
                        void *temp;
                        size += size >> 1;
                        temp = realloc(array, size);
                        if (temp == NULL)
                            goto cleanup;
                        array = temp;
                    }
                    array[ptr++] = mask == 255? data: alphabet[data];
                }
            }
        }

        result = array;
        array = NULL;
        result[ptr] = 0;
    }

cleanup:
    free(array);
    bi_destroy(cpy);
    return bis_new(result, ptr);
}

/* prints bi to string, with specified base
 * returns NULL on out of memory condition
 * releases ownership of returned string */
bigint_string bi_sprint(const bigint *bi, size_t base)
{
    char *result = NULL;
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    bigint *cpy = NULL;
    size_t size = 100, ptr = 0;
    bi_signed_leaf r;
    char *array;

    array = malloc(size);
    if (array == NULL) goto cleanup;

    if (bi->sign)
        array[ptr++] = '-';

    if (bi_is_zero(bi) || base < 2)
    {
        result = array;
        if (base >= 2)
            result[ptr++] = '0';
        result[ptr] = 0;
        array = NULL;
    }
    else if (ispow2(base))
        return bi_sprint_basepow2(bi, array, ptr, size, lg2(base));
    else
    {
        cpy = bi_copy(bi);
        if (cpy == NULL) goto cleanup;
        cpy->sign = 0;

        while (!bi_is_zero(cpy))
        {
            if (ptr == size)
            {
                void *temp;
                size += size >> 1;
                temp = realloc(array, size);
                if (temp == NULL)
                    goto cleanup;
                array = temp;
            }
            if ((cpy = bi_divmod_immediate_assign(cpy, base, &r)) == NULL) goto cleanup;
            array[ptr++] = r <= 36? alphabet[r]: r;
        }

        result = malloc(size+1);
        if (result == NULL)
            goto cleanup;

        size = 0;
        while (ptr > 0)
            result[size++] = array[--ptr];
        result[size] = 0;
        ptr = size;
    }

cleanup:
    free(array);
    bi_destroy(cpy);
    return bis_new(result, ptr);
}

void bi_print(const bigint *bi, size_t base)
{
    bi_fprint(stdout, bi, base);
}

/* scans bi from string str, with specified base */
/* returns -1 and destroys bi if out of memory */
/* returns number of digits scanned otherwise */
int bi_sscan(const char *str, bigint *bi, size_t base)
{
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    const char *p;
    int c = ' ', hasdigits = 0, neg = 0;

    if (base <= 36)
    {
        while (isspace(c))
            c = *str++;
        if (!c) return hasdigits;

        bi_clear(bi);
        if (c == '-' || c == '+')
        {
            neg = c == '-';
            c = *str++;
        }

        if (base < 2) return 0;
    }
    else
    {
        bi_clear(bi);
        c = *str++;
    }

    while (1)
    {
        if (!c) break;

        if (base > 36 || (p = strchr(alphabet, tolower(c))) != NULL)
        {
            ++hasdigits;
            if (bi_mul_immediate_assign(bi, base) == NULL) return -1;
            if (base > 36)
            {
                if (bi_add_immediate_assign(bi, (int) (c % base)) == NULL)
                    return -1;
            }
            else
            {
                c = p-alphabet;
                if ((size_t) c > base)
                    break;
                if (bi_add_immediate_assign(bi, c) == NULL)
                    return -1;
            }
        }
        else
            break;
        c = *str++;
    }

    bi->sign = neg && !bi_is_zero(bi);

    return hasdigits;
}

/* scans bi from string str, with specified base and string length */
/* returns -1 and destroys bi if out of memory */
/* returns number of digits scanned otherwise */
int bi_sscan_n(const char *str, size_t len, bigint *bi, size_t base)
{
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    const char *p;
    int c = ' ', hasdigits = 0, neg = 0;

    if (len == 0) return 0;

    if (base <= 36)
    {
        while (isspace(c) && len > 0)
        {
            c = *str++;
            --len;
        }
        ++len;
        if (!c) return hasdigits;

        bi_clear(bi);
        if (len > 0 && (c == '-' || c == '+'))
        {
            neg = c == '-';
            c = *str++;
            --len;
        }

        if (base < 2 || len == 0) return 0;
    }
    else
    {
        bi_clear(bi);
        c = *str++;
    }

    while (len-- > 0)
    {
        if (base > 36 || (p = strchr(alphabet, tolower(c))) != NULL)
        {
            ++hasdigits;
            if (bi_mul_immediate_assign(bi, base) == NULL) return -1;
            if (base > 36)
            {
                if (bi_add_immediate_assign(bi, (int) (c % base)) == NULL)
                    return -1;
            }
            else
            {
                c = p-alphabet;
                if ((size_t) c > base)
                    break;
                if (bi_add_immediate_assign(bi, c) == NULL)
                    return -1;
            }
        }
        else
            break;
        c = *str++;
    }

    bi->sign = neg && !bi_is_zero(bi);

    return hasdigits;
}

/* scans bi from stream f, with specified base */
/* returns -1 and destroys bi if out of memory */
int bi_fscan(FILE *f, bigint *bi, size_t base)
{
    const char alphabet[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    const char *p;
    int c = ' ', hasdigits = 0, neg = 0;

    if (base <= 36)
    {
        while (isspace(c))
            c = fgetc(f);
        if (c == EOF) return hasdigits;

        bi_clear(bi);
        if (c == '-' || c == '+')
        {
            neg = c == '-';
            c = fgetc(f);
        }

        if (base < 2) return 0;
    }
    else
        bi_clear(bi);

    while (1)
    {
        if (c == EOF) break;
        if (base > 36 || (p = strchr(alphabet, tolower(c))) != NULL)
        {
            hasdigits = 1;
            if (bi_mul_immediate_assign(bi, base) == NULL) return -1;
            if (base > 36)
            {
                if (bi_add_immediate_assign(bi, (int) (c % base)) == NULL)
                    return -1;
            }
            else
            {
                c = p-alphabet;
                if ((size_t) c > base)
                    break;
                if (bi_add_immediate_assign(bi, c) == NULL)
                    return -1;
            }
        }
        else
        {
            ungetc(c, f);
            break;
        }
        c = fgetc(f);
    }

    bi->sign = neg && !bi_is_zero(bi);

    return hasdigits;
}

/* scans bi from stdin, with specified base */
/* returns -1 and destroys bi if out of memory */
int bi_scan(bigint *bi, size_t base)
{
    return bi_fscan(stdin, bi, base);
}

/* destroys a previously allocated bigint */
void bi_destroy(bigint *bi)
{
    if (bi != NULL)
        free(bi->data);
    free(bi);
}

/* initializes a bigint string */
bigint_string bis_new(char *s, size_t len)
{
    bigint_string r;
    r.string = s;
    r.len = len;
    return r;
}

/* destroys a previously allocated bigint string */
void bis_destroy(bigint_string s)
{
    free(s.string);
}
