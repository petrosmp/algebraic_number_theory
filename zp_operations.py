def mulZp_scalar(x: int, y: int, p: int) -> int:
    """Computes the scalar multiplication x*y mod p"""
    return (x * y) % p


def sumZp_list(x: list[int], p: int) -> int:
    """Computes the sum of the elements in the given list mod p"""
    return sum(x) % p


def sumZp(x: list[int], y: list[int], p: int) -> list[int]:
    """
    Calculates x+y mod p. Assumes x,y are of the same length, returns list of the same length.
    Also assumes that all elements of x, y are in Zp.
    """
    x, y = prokroustis(x, y)
    res = []
    for i in range(len(x)):
        res.append((x[i] + y[i]) % p)
    return res


def difZp(x: list[int], y: list[int], p: int) -> list[int]:
    """
    Calculates x-y mod p. Assumes x,y are of the same length, returns list of the same length.
    Also assumes that all elements of x, y are in Zp.
    """
    x, y = prokroustis(x, y)
    res = []
    for i in range(len(x)):
        res.append((x[i] - y[i]) % p)
    return res


def oppZp(x: list[int], p: int) -> list[int]:
    """Calculates -x mod p. Assumes all elements of x are in Zp"""
    return [-i % p for i in x]


def mulZp(x: list[int], y: int, p: int) -> list[int]:
    """Calculates x*y mod p. Assumes y is in Zp, as are all elements of x."""
    return [(i * y) % p for i in x]


def convZp(x: list[int], y: list[int], p: int):
    """Calculates the convolution of 2 polynomials in Zp[x]. Assumes x, y are of the same length."""
    res = []
    x, y = prokroustis(x, y)
    y = [0 for _ in y] + y[::-1]
    for i in range(len(y) - 1):
        x_part = x[: i + 1]  # up to the i'th element (i=0 -> only the first, i>=len -> all)
        y_part = y[-(i + 1) :]  # from the i'th last element onwards (i=0->only the last, i=len -> all)
        res.append(sumZp_list([mulZp_scalar(k, l, p) for k, l in zip(x_part, y_part)], p))
    return res


def invZp(x: int, p: int) -> int | None:
    """
    Calculates the inverse of x in Zp. Returns `None` if it does not exist, i.e. if
    `gcd(x, p) != 1`.

    Basically solves the following congruence: `ax = 1 mod p`, `a = ?` by using the
    extended euclidean algorithm to express gcd(x, p) as a linear combination of x and p.
      - Note that c_b is the answer (corresponds to c_x, since x is less than p and when
        calculating gcd(a,b) b is always the smallest).
    """
    gcd, _, c_x = extended_euclid(x, p)
    if gcd != 1:
        return None
    return c_x % p


def divZp(x, y: int, p) -> list[int] | None:
    """
    Calculates x/y mod Zp. Division is just multiplication with the inverse, so this is what this does

    If the inverse of y does not exist, returns None.
    """
    inv_y = invZp(y, p)
    if not inv_y:
        return None
    return mulZp(x, inv_y, p)


def extended_euclid(a: int, b: int) -> tuple[int, int, int]:
    """
    Calculates the GCD of the given integers and the coefficients of the linear
    combination of them that is equal to it.

    Returns a tuple `(gcd, c_a, c_b)`, with `c_a`, `c_b` such that
    `a*c_a + b*c_b = gcd`.

    The algorithm implemented here is derived in the book as "Algorithm 4".
    """

    a, b = max(a, b), min(a, b)

    r_prev, r = a, b
    x_prev, x = 1, 0
    y_prev, y = 0, 1

    while r != 0:
        q = r_prev // r
        r_prev, r = r, r_prev - q * r
        x_prev, x = x, x_prev - q * x
        y_prev, y = y, y_prev - q * y

    return r_prev, x_prev, y_prev


def polynomial_division_Zp(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    """
    Performs the division a(x) / b(x), where `a` and `b` are polynomials in the integral domain
    `Z_p[x]` (`p` is assumed to be prime).

    a and b are represented by their coefficients (they are basically vectors in `Z^deg(a)` and
    `Z^deg(b)` respectively).

    `a` is assumed to be of greater degree than `b`, and the index of an element specifies the power
    of x of which the element is the coefficient. Trailing 0's are ignored.
    For example, `[4, 6, 0, 1, 3, 5, 0]` represents the polynomial `f(x) = 4 + 6x + x^3 + 3x^4 + 5x^5`.

    Returns a tuple:
        `(q, r)`
    """
    a = remove_trailing_zeros(a)
    b = remove_trailing_zeros(b)

    deg_b = len(b)
    b_leading_coeff = b[-1]

    q = [0 for _ in range(len(a) + deg_b)]  # start it from all 0 so then we can add instead of appending
    r = []

    while deg_b <= len(a):
        # find k: k*b_leading_coeff = a_leading_coeff (mod p) so we can cancel it out
        a_leading_coeff = a[-1]
        k = mulZp_scalar(invZp(b_leading_coeff, p), a_leading_coeff, p)

        # create a polynomial with leading coefficient k and the appropriate degree (i.e. deg(a) - deg(b))
        deg_k = len(a) - deg_b
        shifted_k = [0] * deg_k + [k]
        q = sumZp(q, pad_with_zeros(shifted_k, len(q)), p)

        # multiply b by the quotient we found
        mutliplied = [0]*deg_k + mulZp(b, k, p) # this could also be convZp(b, shifted_k, p)

        # subtract the product from the dividend
        r = difZp(a, mutliplied, p)

        # continue with the remainder as the new dividend
        a = remove_trailing_zeros(r)

    return remove_trailing_zeros(q), r


def remove_trailing_zeros(x: list[int]) -> list[int]:
    while x and x[-1] == 0:
        x.pop()
    return x


def pad_with_zeros(x: list[int], desired_length: int) -> list[int]:
    """
    Append as many zeros as needed to the given list so that it reaches the desired length.
    """
    num_zeros = desired_length - len(x)
    return x + [0 for _ in range(num_zeros)]


def ext_euc_alg_poly(
    a: list[int], b: list[int], p: int, verbose: bool = False
) -> tuple[list[str], list[str], list[str]]:
    """
    Calculates the GCD of the given polynomials and the coefficients of the linear combination of them that
    is equal to it.

    Returns a tuple `(gcd, non_monic, c_a, c_b)`, with `c_a`, `c_b` such that `a*c_a + b*c_b = gcd`.
    `non_monic` is there for debugging purposes, it is the non-monic associate of `gcd` that resulted from the
    algorithm and was then made monic. If the result of the algorithm is already monic, `non_monic` is `None`.

    If the `verbose` flag is set, a summary of each step of the algorithm is printed.

    The algorithm is the same as the one for integers, only the operations are now done on polynomials.
    """
    a, b = max(a, b), min(a, b)

    r_prev, r = remove_trailing_zeros(a), remove_trailing_zeros(b)
    x_prev, x = [1], [0]
    y_prev, y = [0], [1]

    i = 1

    while r:
        q, _r = polynomial_division_Zp(r_prev, r, p)
        r_prev, r = r, _r
        x_prev, x = x, difZp(x_prev, convZp(q, x, p), p)
        y_prev, y = y, difZp(y_prev, convZp(q, y, p), p)

        if verbose:
            print(f"step {i}: {to_polynomial(r)} = {to_polynomial(x)}*a + {to_polynomial(y)}*b")
            print(f"\tto verify, x*a = {to_polynomial(remove_trailing_zeros(convZp(x, a, p)))}")
            print(f"\tto verify, y*b = {to_polynomial(remove_trailing_zeros(convZp(y, b, p)))}")
            print(
                f"\tto verify, sum = {to_polynomial(remove_trailing_zeros(sumZp(convZp(x, a, p), convZp(y, b, p), p)))}"
            )
            print(f"\tto verify,  r  = {to_polynomial(remove_trailing_zeros(r))}")
            print()
            i += 1

    # check if the result is monic, and if not make it
    result_leading_coeff = r_prev[-1]
    non_monic = None
    if result_leading_coeff != 1:
        non_monic = r_prev
        inv = invZp(result_leading_coeff, p)
        r_prev = mulZp(r_prev, inv, p)
        x_prev = mulZp(x_prev, inv, p)
        y_prev = mulZp(y_prev, inv, p)

    return r_prev, non_monic, x_prev, y_prev


def to_polynomial(x: list[int]) -> str:
    if not x:
        return "0"

    res = ""
    first = True
    for i, coeff in enumerate(x):
        if coeff != 0:
            if i == 0:
                res += f"{coeff}"
                first = False
            elif i == 1:
                if first:
                    first = False
                else:
                    res += " + "
                res += f"{coeff if coeff != 1 else ''}x"
            else:
                if first:
                    first = False
                else:
                    res += " + "
                res += f"{coeff if coeff != 1 else ''}x^{i}"
    return res


def prokroustis(x: list[int], y: list[int]) -> tuple[list[int], list[int]]:
    """
    Makes the given arrays the same length by appending zeros to the smallest one.
    """
    l_x = len(x)
    l_y = len(y)
    if l_x > l_y:
        y = pad_with_zeros(y, l_x)
    else:
        x = pad_with_zeros(x, l_y)

    return x, y
