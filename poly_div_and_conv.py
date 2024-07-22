from zp_operations import *

a = [1, 2, 0, 2, 1]
b = [1, 0, 2, 1]
p = 3

print(f"a(x) = {to_polynomial(a)}")
print(f"b(x) = {to_polynomial(b)}")

q, r = polynomial_division_Zp(a, b, p)
print(f"a(x) = ({to_polynomial(q)}) * b(x) + {to_polynomial(r)}")

gcd, non_monic, c_a, c_b = ext_euc_alg_poly(a, b, p)
print(f"gcd = {to_polynomial(gcd)} (the monic associate of {to_polynomial(non_monic)})")
print(f"a(x) coefficient: {to_polynomial(c_a)}")
print(f"b(x) coefficient: {to_polynomial(c_b)}")

print("\nVerification:")
print(f"a(x) * c_a(x) = {to_polynomial(convZp(a, c_a, p))}")
print(f"b(x) * c_b(x) = {to_polynomial(convZp(b, c_b, p))}")
print(f"b(x)*c_b(x) + a*c_a(x) = {to_polynomial(sumZp(convZp(a, c_a, p), convZp(b, c_b, p), p))}")
