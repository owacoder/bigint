# bigint
Arbitrary precision math library

To use this library, simply include the relevant source and header files in your project.
The bigfrac library depends on the bigint library, but the bigint library only depends on
the standard C (or C++, if compiled with __cplusplus defined) library.

If using the raw C interface, return values should be checked to verify that the library
has not run out of memory. The C++ interface will throw an out-of-memory error when necessary.

A small example program using the C interface:

    int main()
    {
        bigint *b = bi_new();
        size_t i;
        
        if (b == NULL || bi_assign(b, 2) == NULL)
            return 1;
        
        for (i = 0; i < 10; ++i)
        {
            if (bi_square_assign(b) == NULL)
                return 1;
            bi_print(b, 10);
            puts("");
        }
        
        bi_destroy(b);
        return 0;
    }
    
Here is the same example program written with the C++ interface:

    int main()
    {
        try {
            Bigint b = 2;

            for (size_t i = 0; i < 10; ++i)
            {
                b.square();
                std::cout << b.toString() << std::endl;
            }
        } catch (Bigint::out_of_memory) {
            return 1;
        }
        
        return 0;
    }
    
# Features

 - Basic math: Addition, subtraction, multiplication, division, modulus.
 - Bitwise math: AND, OR, XOR, SHL, SHR (logical).
 - Comparison: By value, by magnitude.
 - Exponential functions: square (faster than basic multiplication), square root, raise to power, log2.
 - Optional exponential functions (without BIGINT_DISABLE_LIBMATH defined): log10, logn (any integral base n >= 2)
 - Additional functions: abs, negate, factorial, GCD.
 - Input/Output: scan in base 2-256 from `FILE*` or `char*`, print in base 2-256 to `FILE*` or `char*`.
    
# Compile flags

By default, the minimum number of leaves is 8. This can be changed externally by pre-defining `BIGINT_MIN_LEAFS` to a different number. Note this value must be greater than zero, or the library will not build.

The leaf size should default to the native word size of the target platform. However, this size can be overridden by pre-defining `BIGINT_WORD_SIZE` to 8, 16, or 32. On GCC (or MinGW) target platforms, `BIGINT_WORD_SIZE` can also be set to 64.

The thresholds for Karatsuba and Toom-Cook multiplication can also be modified. They are set to reasonable (though arbitrary) defaults, and can be set by pre-defining `BIGINT_KARATSUBA_MIN_LEAFS` and `BIGINT_TOOM_COOK_MIN_LEAFS`, respectively.

On a GCC (or MinGW) target, pthreads are enabled by default. They can be turned off by pre-defining `BIGINT_DISABLE_PTHREADS`.

If on a target without floating-point math support, define `BIGINT_DISABLE_LIBMATH`. This will prevent the build from using any floating-point math, and will remove (i.e. not build) the `log10`, `log10_approx`, `logn`, and `logn_approx` functions.
