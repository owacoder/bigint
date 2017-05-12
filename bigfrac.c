#include "bigfrac.h"
#include "general.h"

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>

/* creates a new bigfrac initialized to 0/1
 * returns NULL on out of memory */
bigfrac *bf_new()
{
    bigfrac *ptr = malloc(sizeof(bigfrac));
    if (ptr != NULL)
    {
        ptr->n = bi_new();
        ptr->d = bi_new();
        if (ptr->n == NULL || ptr->d == NULL || bi_assignu(ptr->d, 1) == NULL)
            goto cleanup;
    }
    return ptr;

cleanup:
    if (ptr != NULL)
    {
        bi_destroy(ptr->n);
        bi_destroy(ptr->d);
        free(ptr);
    }
    return NULL;
}

bigfrac *bf_new_value(bi_signed_leaf value)
{
    return bf_new_frac_init(bi_new_value(value), bi_new_valueu(1));
}

bigfrac *bf_new_valueu(bi_leaf value)
{
    return bf_new_frac_init(bi_new_valueu(value), bi_new_valueu(1));
}

bigfrac *bf_new_valuel(bi_intmax value)
{
    return bf_new_frac_init(bi_new_valuel(value), bi_new_valueu(1));
}

bigfrac *bf_new_valuelu(bi_uintmax value)
{
    return bf_new_frac_init(bi_new_valuelu(value), bi_new_valueu(1));
}

bigfrac *bf_new_frac_value(bi_signed_leaf numer, bi_signed_leaf denom)
{
    return bf_new_frac_init(bi_new_value(numer), bi_new_value(denom));
}

bigfrac *bf_new_frac_valueu(bi_leaf numer, bi_leaf denom)
{
    return bf_new_frac_init(bi_new_valueu(numer), bi_new_valueu(denom));
}

bigfrac *bf_new_frac_valuel(bi_intmax numer, bi_intmax denom)
{
    return bf_new_frac_init(bi_new_valuel(numer), bi_new_valuel(denom));
}

bigfrac *bf_new_frac_valuelu(bi_uintmax numer, bi_uintmax denom)
{
    return bf_new_frac_init(bi_new_valuelu(numer), bi_new_valueu(denom));
}

/* creates a new bigfrac initialized to specified values
 * returns NULL on out of memory */
bigfrac *bf_new_frac(const bigint *numer, const bigint *denom)
{
    bigfrac *ptr = malloc(sizeof(bigfrac));
    if (ptr != NULL)
    {
        ptr->n = bi_copy(numer);
        ptr->d = bi_copy(denom);
        if (ptr->n == NULL || ptr->d == NULL)
            goto cleanup;
    }
    return ptr;

cleanup:
    if (ptr != NULL)
    {
        bi_destroy(ptr->n);
        bi_destroy(ptr->d);
        free(ptr);
    }
    return NULL;
}

/* creates a new bigfrac initialized to the specified values, taking ownership
 * returns NULL and destroys parameters on out of memory */
bigfrac *bf_new_frac_init(bigint *numer, bigint *denom)
{
    bigfrac *ptr;

    if (numer == NULL || denom == NULL)
    {
        bi_destroy(numer);
        bi_destroy(denom);
        return NULL;
    }

    ptr = malloc(sizeof(bigfrac));
    if (ptr != NULL)
    {
        ptr->n = numer;
        ptr->d = denom;
    }
    else
    {
        bi_destroy(numer);
        bi_destroy(denom);
    }
    return ptr;
}

/* copies fraction
 * returns NULL on out of memory */
bigfrac *bf_copy(const bigfrac *bf)
{
    bigfrac *ptr = malloc(sizeof(bigfrac));
    if (ptr != NULL)
    {
        ptr->n = bi_copy(bf->n);
        ptr->d = bi_copy(bf->d);
        if (ptr->n == NULL || ptr->d == NULL)
            goto cleanup;
    }
    return ptr;

cleanup:
    if (ptr != NULL)
    {
        bi_destroy(ptr->n);
        bi_destroy(ptr->d);
        free(ptr);
    }
    return NULL;
}

/* copies src to dst
 * returns NULL and destroys dst on out of memory */
bigfrac *bf_copy_to(bigfrac *dst, const bigfrac *src)
{
    if (dst == src)
        return dst;

    dst->n = bi_copy(src->n);
    dst->d = bi_copy(src->d);
    if (dst->n == NULL || dst->d == NULL)
    {
        bf_destroy(dst);
        return NULL;
    }
    return dst;
}

/* copies numer and denom to dst
 * returns NULL and destroys dst on out of memory */
bigfrac *bf_copy_frac_to(bigfrac *dst, const bigint *numer, const bigint *denom)
{
    dst->n = bi_copy_to(dst->n, numer);
    dst->d = bi_copy_to(dst->d, denom);
    if (dst->n == NULL || dst->d == NULL)
    {
        bf_destroy(dst);
        return NULL;
    }
    return dst;
}

/* copies numer and denom to dst, taking ownership of numer and denom
 * returns NULL and destroys dst and parameters on out of memory */
bigfrac *bf_copy_frac_init_to(bigfrac *dst, bigint *numer, bigint *denom)
{
    if (numer == NULL || denom == NULL)
    {
        bi_destroy(numer); dst->n=NULL;
        bi_destroy(denom); dst->d=NULL;
        bf_destroy(dst);
        return NULL;
    }
    if (dst->n != numer) {bi_destroy(dst->n); dst->n = numer;}
    if (dst->d != denom) {bi_destroy(dst->d); dst->d = denom;}
    return dst;
}

/* copies numer to bf
 * returns NULL and destroys bf on out of memory */
bigfrac *bf_copy_numer_to(bigfrac *bf, const bigint *numer)
{
    bf->n = bi_copy_to(bf->n, numer);
    if (bf->n == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_copy_denom_to(bigfrac *bf, const bigint *denom)
{
    bf->d = bi_copy_to(bf->d, denom);
    if (bf->d == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_clear(bigfrac *bf)
{
    bi_clear(bf->n);
    bi_clear(bf->d);
    return bf;
}

bigfrac *bf_assign(bigfrac *bf, bi_signed_leaf value)
{
    if ((bf->n = bi_assign(bf->n, value)) == NULL ||
        (bf->d = bi_assignu(bf->d, 1)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_assignu(bigfrac *bf, bi_leaf value)
{
    if ((bf->n = bi_assignu(bf->n, value)) == NULL ||
        (bf->d = bi_assignu(bf->d, 1)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_assignl(bigfrac *bf, bi_intmax value)
{
    if ((bf->n = bi_assignl(bf->n, value)) == NULL ||
        (bf->d = bi_assignu(bf->d, 1)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_assignlu(bigfrac *bf, bi_uintmax value)
{
    if ((bf->n = bi_assignlu(bf->n, value)) == NULL ||
        (bf->d = bi_assignu(bf->d, 1)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_frac_assign(bigfrac *bf, bi_signed_leaf numer, bi_signed_leaf denom)
{
    if ((bf->n = bi_assign(bf->n, numer)) == NULL ||
        (bf->d = bi_assign(bf->d, denom)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_frac_assignu(bigfrac *bf, bi_leaf numer, bi_leaf denom)
{
    if ((bf->n = bi_assignu(bf->n, numer)) == NULL ||
        (bf->d = bi_assignu(bf->d, denom)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_frac_assignl(bigfrac *bf, bi_intmax numer, bi_intmax denom)
{
    if ((bf->n = bi_assignl(bf->n, numer)) == NULL ||
        (bf->d = bi_assignl(bf->d, denom)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

bigfrac *bf_frac_assignlu(bigfrac *bf, bi_uintmax numer, bi_uintmax denom)
{
    if ((bf->n = bi_assignlu(bf->n, numer)) == NULL ||
        (bf->d = bi_assignlu(bf->d, denom)) == NULL)
    {
        bf_destroy(bf);
        return NULL;
    }
    return bf;
}

/* returns 0 if bf == bf2, 1 if bf != bf2, 2 if out of memory */
int bf_cmp_eq(const bigfrac *bf, const bigfrac *bf2)
{
    bigint *a = bi_mul(bf->n, bf2->d);
    bigint *b = bi_mul(bf->d, bf2->n);
    int r;

    if (a == NULL || b == NULL)
    {
        bi_destroy(a);
        bi_destroy(b);
        return 2;
    }

    r = bi_cmp(a, b);
    bi_destroy(a);
    bi_destroy(b);

    return r != 0;
}

/* returns -1 if bf < bf2, 0 if bf == bf2, 1 if bf > bf2, 2 if out of memory */
int bf_cmp(const bigfrac *bf, const bigfrac *bf2)
{
    bigint *a = bi_mul(bf->n, bf2->d);
    bigint *b = bi_mul(bf->d, bf2->n);
    int r;

    if (a == NULL || b == NULL)
    {
        bi_destroy(a);
        bi_destroy(b);
        return 2;
    }

    a->sign = bi_is_negative(bf->d)? !bf->n->sign: bf->n->sign;
    b->sign = bi_is_negative(bf2->d)? !bf2->n->sign: bf2->n->sign;

    r = bi_cmp(a, b);
    bi_destroy(a);
    bi_destroy(b);

    return r;
}

/* returns -1 if bf < 0, 0 if bf == 0, 1 if bf > 0 */
int bf_cmp_zero(const bigfrac *bf)
{
    if (bf_is_negative(bf))
        return -1;
    else if (bf_is_zero(bf))
        return 0;
    else
        return 1;
}

int bf_is_negative(const bigfrac *bf)
{
    return bi_is_negative(bf->n) != bi_is_negative(bf->d);
}

int bf_is_one(const bigfrac *bf)
{
    return bi_cmp(bf->n, bf->d) == 0 && !bi_is_zero(bf->d);
}

int bf_is_zero(const bigfrac *bf)
{
    return bi_is_zero(bf->n) && !bi_is_zero(bf->d);
}

int bf_is_undefined(const bigfrac *bf)
{
    return bi_is_zero(bf->d);
}

int bf_is_power_of_2(const bigfrac *bf)
{
    return bi_is_power_of_2(bf->n) && bi_is_power_of_2(bf->d);
}

bigfrac *bf_negate(const bigfrac *bf)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_negate_assign(result);
}

bigfrac *bf_negate_assign(bigfrac *bf)
{
    bi_negate_assign(bf->n);
    return bf;
}

bigfrac *bf_abs(const bigfrac *bf)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_abs_assign(result);
}

bigfrac *bf_abs_assign(bigfrac *bf)
{
    bi_abs_assign(bf->n);
    bi_abs_assign(bf->d);
    return bf;
}

bigfrac *bf_reciprocal(const bigfrac *bf)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_reciprocal_assign(result);
}

bigfrac *bf_reciprocal_assign(bigfrac *bf)
{
    bigint *temp = bf->n;
    bf->n = bf->d;
    bf->d = temp;
    return bf;
}

bigfrac *bf_add(const bigfrac *bf, const bigfrac *bf2)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_add_assign(result, bf2);
}

bigfrac *bf_add_assign(bigfrac *bf, const bigfrac *bf2)
{
    bigint *n1 = bi_mul(bf->n, bf2->d);
    bigint *n2 = bi_mul(bf->d, bf2->n);
    bigint *n = NULL;
    if (n1 == NULL || n2 == NULL)
        goto cleanup;

    n = bi_add(n1, n2);
    if (n == NULL)
        goto cleanup;
    bi_destroy(n1); n1 = NULL;
    bi_destroy(n2); n2 = NULL;
    return bf_copy_frac_init_to(bf, n, bi_mul(bf->d, bf2->d));

cleanup:
    bi_destroy(n1);
    bi_destroy(n2);
    bi_destroy(n);
    bf_destroy(bf);
    return NULL;
}

bigfrac *bf_add_immediate(const bigfrac *bf, bi_signed_leaf val)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_add_immediate_assign(result, val);
}

bigfrac *bf_add_immediate_assign(bigfrac *bf, bi_signed_leaf val)
{
    if (bi_add_immediate_assign(bf->n, val) == NULL)
    {
        bi_destroy(bf->d);
        bf->n = bf->d = NULL;
        bf_destroy(bf);
    }
    return bf;
}

bigfrac *bf_sub(const bigfrac *bf, const bigfrac *bf2)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_sub_assign(result, bf2);
}

bigfrac *bf_sub_assign(bigfrac *bf, const bigfrac *bf2)
{
    bigint *n1 = bi_mul(bf->n, bf2->d);
    bigint *n2 = bi_mul(bf->d, bf2->n);
    bigint *n = NULL;
    if (n1 == NULL || n2 == NULL)
        goto cleanup;

    n = bi_sub(n1, n2);
    if (n == NULL)
        goto cleanup;
    bi_destroy(n1); n1 = NULL;
    bi_destroy(n2); n2 = NULL;
    return bf_copy_frac_init_to(bf, n, bi_mul(bf->d, bf2->d));

cleanup:
    bi_destroy(n1);
    bi_destroy(n2);
    bi_destroy(n);
    return NULL;
}

bigfrac *bf_sub_immediate(const bigfrac *bf, bi_signed_leaf val)
{
    bigfrac *result = bf_copy(bf);
    if (result == NULL) return NULL;
    return bf_sub_immediate_assign(result, val);
}

bigfrac *bf_sub_immediate_assign(bigfrac *bf, bi_signed_leaf val)
{
    if (bi_sub_immediate_assign(bf->n, val) == NULL)
    {
        bi_destroy(bf->d);
        bf->n = bf->d = NULL;
        bf_destroy(bf);
    }
    return bf;
}

bigfrac *bf_mul(const bigfrac *bf, const bigfrac *bf2)
{
    return bf_new_frac_init(bi_mul(bf->n, bf2->n), bi_mul(bf->d, bf2->d));
}

bigfrac *bf_mul_assign(bigfrac *bf, const bigfrac *bf2)
{
    return bf_copy_frac_init_to(bf, bi_mul(bf->n, bf2->n), bi_mul(bf->d, bf2->d));
}

bigfrac *bf_mul_immediateu(const bigfrac *bf, bi_leaf val)
{
    return bf_new_frac_init(bi_mul_immediateu(bf->n, val), bi_copy(bf->d));
}

bigfrac *bf_mul_immediateu_assign(bigfrac *bf, bi_leaf val)
{
    return bf_copy_frac_init_to(bf, bi_mul_immediateu_assign(bf->n, val), bf->d);
}

bigfrac *bf_mul_immediate(const bigfrac *bf, bi_signed_leaf val)
{
    return bf_new_frac_init(bi_mul_immediate(bf->n, val), bi_copy(bf->d));
}

bigfrac *bf_mul_immediate_assign(bigfrac *bf, bi_signed_leaf val)
{
    return bf_copy_frac_init_to(bf, bi_mul_immediate_assign(bf->n, val), bf->d);
}

bigfrac *bf_square(const bigfrac *bf)
{
    return bf_new_frac_init(bi_square(bf->n), bi_square(bf->d));
}

bigfrac *bf_square_assign(bigfrac *bf)
{
    return bf_copy_frac_init_to(bf, bi_square_assign(bf->n), bi_square_assign(bf->d));
}

bigfrac *bf_sqrt(const bigfrac *bf)
{
    return bf_new_frac_init(bi_sqrt(bf->n), bi_sqrt(bf->d));
}

bigfrac *bf_sqrt_assign(bigfrac *bf)
{
    if ((bf = bf_reduce(bf)) == NULL)
        return NULL;
    return bf_copy_frac_init_to(bf, bi_sqrt_assign(bf->n), bi_sqrt_assign(bf->d));
}

bigfrac *bf_div(const bigfrac *bf, const bigfrac *bf2)
{
    return bf_new_frac_init(bi_mul(bf->n, bf2->d), bi_mul(bf->d, bf2->n));
}

bigfrac *bf_div_assign(bigfrac *bf, const bigfrac *bf2)
{
    return bf_copy_frac_init_to(bf, bi_mul(bf->n, bf2->d), bi_mul(bf->d, bf2->n));
}

bigfrac *bf_div_immediateu(const bigfrac *bf, bi_leaf val)
{
    return bf_new_frac_init(bi_copy(bf->n), bi_mul_immediateu(bf->d, val));
}

bigfrac *bf_div_immediateu_assign(bigfrac *bf, bi_leaf val)
{
    return bf_copy_frac_init_to(bf, bf->n, bi_mul_immediateu(bf->d, val));
}

bigfrac *bf_div_immediate(const bigfrac *bf, bi_signed_leaf val)
{
    return bf_new_frac_init(bi_copy(bf->n), bi_mul_immediate(bf->d, val));
}

bigfrac *bf_div_immediate_assign(bigfrac *bf, bi_signed_leaf val)
{
    return bf_copy_frac_init_to(bf, bf->n, bi_mul_immediate(bf->d, val));
}

bigfrac *bf_reduce(bigfrac *bf)
{
    bigint *gcd;
    size_t negative = bf->n->sign != bf->d->sign;

    if (bi_is_zero(bf->n) || bi_is_zero(bf->d))
        return bf;
    else if (bi_cmp_mag(bf->n, bf->d) == 0)
    {
        if ((bf->n = bi_assign(bf->n, negative? -1: 1)) == NULL ||
            (bf->d = bi_assignu(bf->d, 1)) == NULL)
            goto cleanup;
        return bf;
    }

    gcd = bi_gcd(bf->n, bf->d);
    if (gcd == NULL) goto cleanup;

    bf->n = bi_div_assign(bf->n, gcd);
    bf->d = bi_div_assign(bf->d, gcd);
    bi_destroy(gcd); gcd = NULL;
    if (bf->n != NULL && bf->d != NULL)
    {
        // fix sign
        bf->n->sign = negative;
        bf->d->sign = 0;
        return bf;
    }

cleanup:
    bi_destroy(gcd);
    bf_destroy(bf);
    return NULL;
}

int bf_sscan(const char *s, bigfrac *bf, size_t base)
{
    int n = bi_sscan(s, bf->n, base), n2 = 0;

    if (n < 0)
    {
        bf->n = NULL;
        bf_destroy(bf);
        return n;
    }

    if (s[n] == '/')
    {
        s += n + 1;
        n2 = bi_sscan(s, bf->d, base);
        if (n2 < 0)
        {
            bf->d = NULL;
            bf_destroy(bf);
            return n2;
        }
        ++n2; // add length of '/' character
    }
    else if ((bf->d = bi_assign(bf->d, 1)) == NULL)
    {
        bf_destroy(bf);
        return -1;
    }

    return n + n2;
}

int bf_sscan_n(const char *s, size_t len, bigfrac *bf, size_t base)
{
    int n = bi_sscan_n(s, len, bf->n, base), n2 = 0;

    if (n < 0)
    {
        bf->n = NULL;
        bf_destroy(bf);
        return n;
    }
    len -= n;

    if (len > 0 && s[n] == '/')
    {
        s += n + 1;
        n2 = bi_sscan_n(s, --len, bf->d, base);
        if (n2 < 0)
        {
            bf->d = NULL;
            bf_destroy(bf);
            return n2;
        }
        ++n2; // add length of '/' character
    }
    else if ((bf->d = bi_assign(bf->d, 1)) == NULL)
    {
        bf_destroy(bf);
        return -1;
    }

    return n + n2;
}

bigint_string bf_sprint(const bigfrac *bf, size_t base)
{
    bigint_string numer, denom;
    char *temp;

    numer = bi_sprint(bf->n, base);
    denom = bi_sprint(bf->d, base);

    if (numer.string == NULL || denom.string == NULL)
    {
        bis_destroy(numer);
        bis_destroy(denom);
        return bis_new(NULL, 0);
    }

    temp = realloc(numer.string, numer.len + denom.len + 2);
    if (temp == NULL)
    {
        bis_destroy(numer);
        bis_destroy(denom);
        return bis_new(NULL, 0);
    }
    numer.string = temp;

    numer.string[numer.len] = '/';
    memcpy(numer.string + numer.len + 1, denom.string, denom.len);
    numer.string[numer.len + denom.len + 1] = 0;
    bis_destroy(denom);

    return bis_new(numer.string, numer.len + denom.len + 1);
}

void bf_print(const bigfrac *bf, size_t base)
{
    bf_fprint(stdout, bf, base);
}

void bf_fprint(FILE *out, const bigfrac *bf, size_t base)
{
    bi_fprint(out, bf->n, base);
    fputc('/', out);
    bi_fprint(out, bf->d, base);
}

void bf_swap(bigfrac *bfa, bigfrac *bfb)
{
    bigfrac temp = *bfa;
    *bfa = *bfb;
    *bfb = temp;
}

void bf_destroy(bigfrac *bf)
{
    if (bf != NULL)
    {
        bi_destroy(bf->n);
        bi_destroy(bf->d);
    }
    free(bf);
}
