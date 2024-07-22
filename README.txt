                        TEL 511 - Number Theory & Cryptography
                                Spring Semester 2024
                                    Homework 2
                                     Group 6
                                                                         Petros Bimpiris
                                                                    Konstantinos Nikolos

The submission is comprised of the following files:

    - Homework2.pdf: contains the answers to all questions of the exercise statement

    - zp_operations.py: contains the implementation of all the functions described in
      the exercise statement

    - poly_div_conv_gcd.py: contains a demonstration of the following functionalities:
        > polynomial division ("deconvolution")
        > extended euclidean algorithm for polynomials
        > polynomial multiplication (convolution) and addition to verify the result of
          the extended euclidean algorithm

    - construct_galois_field.py: contains a demonstration of the construction of a galois
      field (specifically the one needed for exercise 8 of the problem set) from a given
      prime p and a polynomial f(x) irreducible in Zp[x], using the functions defined
      in zp_operations.py.

      By changing the values of p and f one can easily construct arbitrarily large galois
      fields using the provided functions.

      Note that p is assumed to be prime and f(x) is assumed to be irreducible. If one
      of those assumptions does not hold, no error is raised and the addition/multiplication
      tables are constructed normally, however they do not represent a field, since there
      will be elements that do not have a multiplicative inverse. 

Throughout the code polynomials are represented by an array of coefficients, with the index
of each coefficient representing the power of x to which the coefficient applies. For example,
f(x)  = 1 + x^2 + 3x^4 would be represented as f = [1, 0, 1, 0, 3]. Trailing zeros are
ignored. Coefficient arrays can be transformed to string representations using the
zp_operations.to_polynomial() function, as demonstrated in poly_div_conv_gcd.py.
