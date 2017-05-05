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
bigfrac *bf_new_frac(const bigint *numer, const bigint *denom);
bigfrac *bf_new_frac_init(bigint *numer, bigint *denom);
bigfrac *bf_copy(const bigfrac *bf);
bigfrac *bf_copy_to(bigfrac *dst, const bigfrac *src);
bigfrac *bf_copy_frac_to(bigfrac *dst, const bigint *numer, const bigint *denom);
bigfrac *bf_copy_frac_init_to(bigfrac *dst, bigint *numer, bigint *denom);
bigfrac *bf_copy_numer_to(bigfrac *bf, const bigint *numer);
bigfrac *bf_copy_denom_to(bigfrac *bf, const bigint *denom);
bigfrac *bf_clear(bigfrac *bf);
int bf_cmp_eq(const bigfrac *bf, const bigfrac *bf2);
int bf_cmp(const bigfrac *bf, const bigfrac *bf2);
int bf_is_one(const bigfrac *bf);
int bf_is_zero(const bigfrac *bf);
int bf_is_undefined(const bigfrac *bf);
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
bigfrac *bf_div(const bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_div_assign(bigfrac *bf, const bigfrac *bf2);
bigfrac *bf_div_immediate(const bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_div_immediate_assign(bigfrac *bf, bi_signed_leaf val);
bigfrac *bf_div_immediateu(const bigfrac *bf, bi_leaf val);
bigfrac *bf_div_immediateu_assign(bigfrac *bf, bi_leaf val);
bigfrac *bf_reduce(bigfrac *bf);
bigint_string bf_sprint(const bigfrac *bf, size_t base);
void bf_print(const bigfrac *bf, size_t base);
void bf_fprint(FILE *out, const bigfrac *bf, size_t base);
void bf_swap(bigfrac *bfa, bigfrac *bfb);
void bf_destroy(bigfrac *bf);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // BIGFRAC_H

