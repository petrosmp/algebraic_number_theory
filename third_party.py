"""
A script that uses the galois package (https://pypi.org/project/galois/, version 0.4.1) to
verify the results produced by our own implementation.

Notice that despite being built on top of numpy (which in turn is basically C), the galois
package is noticeably slower than our implementation for some reason.
"""

import galois

gf9 = galois.GF(9, irreducible_poly=[1, 0, 1])
print(gf9.properties, end="\n\n")
print(gf9.arithmetic_table("+"), end="\n\n")
print(gf9.arithmetic_table("*"), end="\n\n")
print(gf9.repr_table(sort="int"))
