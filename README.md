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
            size_t i;

            for (i = 0; i < 10; ++i)
            {
                b.square();
                std::cout << b.toString() << std::endl;
            }
        } catch (Bigint::out_of_memory) {
            return 1;
        }
        
        return 0;
    }
