import math
from sympy import GF, invert, isprime
import numpy as np
from sympy.abc import x
from sympy import ZZ, Poly


def is_2_power(n):
    """Check if a number is a power of 2.

    Determines whether the given integer `n` is a power of 2. This function
    returns `True` for any positive integer that is a power of 2 and `False` otherwise.

    Args:
        n (int): The integer to check.

    Returns:
        bool: `True` if `n` is a power of 2, otherwise `False`.
    """
    return n != 0 and (n & (n - 1) == 0)


def random_poly(length, d, neg_ones_diff=0):
    """Generates a random polynomial with coefficients in Z_2.

    The function creates a polynomial of degree 'length - 1', containing:
    - 'd' positive coefficients (1s),
    - `d + neg_ones_diff` negative coefficients (-1s),
    - Rest are set to Zero (0).
    The coefficients are randomly permuted to ensure varied distribution.

    Args:
        length (int): Total number of coefficients in the polynomial.
                        >= 2 * d + neg_ones_diff.
        d (int): Number of positive coefficients (1s) in the polynomial.
                    >= 0.
        neg_ones_diff (int, optional): Number of negative coefficients (-1s)
                                        beyond 'd' in the polynomial.
                                        >= 0.

    Returns:
        Poly: A polynomial object with coefficients of 0, 1, and -1
                distributed randomly.

    Raises:
        AssertionError: If length is less than 2 * d + neg_ones_diff.
    """
    assert (
        length >= 2 * d + neg_ones_diff
    ), "length should be greater than 2 * d + neg_ones_diff"
    return Poly(
        np.random.permutation(
            np.concatenate(
                (
                    np.zeros(length - 2 * d - neg_ones_diff),
                    np.ones(d),
                    -np.ones(d + neg_ones_diff),
                )
            )
        ),
        x,
    ).set_domain(ZZ)


def invert_poly(f_poly, R_poly, p):
    """Compute the modular inverse of a polynomial over a finite field.

    Calculates the modular inverse of `f_poly` modulo `R_poly` in the field
    defined by modulus `p`. If `p` is prime, the function directly computes
    the inverse in the finite field `GF(p)`. If `p` is a power of 2, the function
    applies Newton's method for refinement to compute the inverse in `GF(2^e)`.

    Args:
        f_poly (Poly): The polynomial to invert.
        R_poly (Poly): The modulus polynomial, which defines the field.
        p (int): The modulus for the finite field. Can be a prime number
                 or a power of 2.

    Returns:
        Poly: The modular inverse of `f_poly` modulo `R_poly`.

    Raises:
        Exception: If `p` is neither prime nor a power of 2, indicating
                    an unsupported modulus for inversion.
    """
    inv_poly = None
    if isprime(p):
        inv_poly = invert(
            f_poly.set_domain(GF(p)), R_poly.set_domain(GF(p)), domain=GF(p)
        )
        return Poly(inv_poly.all_coeffs(), x).set_domain(ZZ)
    elif is_2_power(p):
        inv_poly = invert(
            f_poly.set_domain(GF(2)), R_poly.set_domain(GF(2)), domain=GF(2)
        )
        inv_poly = Poly(inv_poly.all_coeffs(), x).set_domain(ZZ)
        e = int(math.log(p, 2))
        for _ in range(1, e):
            inv_poly = ((2 * inv_poly - f_poly * inv_poly**2) % R_poly).trunc(p)
        return inv_poly

    else:
        raise Exception(f"Cannot invert polynomial in Z_{p}")


def padding_encode(input_arr, block_size):
    """Adds padding to an array to ensure it aligns with a specified block size.

    The function pads `input_arr` with zeros until its length is a multiple of
    `block_size`. If padding is added, an additional block with a sequence of
    ones is appended to denote padding length.

    Args:
        input_arr (np.ndarray): Input array to be padded.
        block_size (int): Desired block size for padding.

    Returns:
        np.ndarray: The padded array with an appended block indicating padding length.
    """
    input_arr = np.asarray(input_arr)
    n = block_size - len(input_arr) % block_size
    if n == block_size:
        return np.pad(input_arr, (0, n), "constant")
    last_block = np.pad(np.ones(n), (block_size - n, 0), "constant")
    return np.concatenate((np.pad(input_arr, (0, n), "constant"), last_block))


def padding_decode(input_arr, block_size):
    """Removes padding from an array encoded with `padding_encode`.

    The function inspects the last block of `input_arr` to identify padding length
    and removes the corresponding number of padded elements.

    Args:
        input_arr (np.ndarray): Padded array from which to remove padding.
        block_size (int): Block size used during padding.

    Returns:
        np.ndarray: Array with padding removed, restoring the original data.
    """
    last_block = input_arr[-block_size:]
    zeros_to_remove = len(np.trim_zeros(last_block))
    return input_arr[: -(block_size + zeros_to_remove)]
