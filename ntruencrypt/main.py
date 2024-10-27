from .cipher import NTRUCipher
from .utils import random_poly, padding_encode, padding_decode
from sympy.abc import x
from sympy import ZZ, Poly
import numpy as np
import math

pub_key = dict()
priv_key = dict()


def generate_keys(N, p, q):
    """Generates and returns public and private keys for NTRU encryption.

    Initializes an NTRU cipher instance with the given parameters, generates
    the public and private keys, and extracts their polynomial coefficients.

    Args:
        N (int): Degree of the polynomial ring modulus.
        p (int): Small modulus, usually a prime.
        q (int): Large modulus, typically much greater than `p`.

    Returns:
        tuple: Contains two dictionaries representing the private and public keys:
            - priv_key (dict): Dictionary with keys 'N', 'p', 'q', 'f', 'f_p',
              representing the private key parameters and polynomial coefficients.
            - pub_key (dict): Dictionary with keys 'N', 'p', 'q', 'h',
              representing the public key parameters and polynomial coefficients.

    Raises:
        Exception: If the random polynomial generation fails after the set attempts.
    """
    ntru = NTRUCipher(N, p, q)
    ntru.generate_random_keys()
    h = np.array(ntru.h_poly.all_coeffs()[::-1])
    f = ntru.f_poly.all_coeffs()[::-1]
    f_p = ntru.f_p_poly.all_coeffs()[::-1]

    priv_key["N"] = N
    priv_key["p"] = p
    priv_key["q"] = q
    priv_key["f"] = ntru.f_coeffs
    priv_key["f_p"] = ntru.f_p_coeffs

    pub_key["N"] = N
    pub_key["p"] = p
    pub_key["q"] = q
    pub_key["h"] = ntru.h_coeffs


def encrypt(msg):
    """Encrypts a message using the NTRU cryptographic scheme.

    The function initializes an NTRU cipher with the public key and encrypts
    the input message. If the message length exceeds `N`, it is padded and split
    into blocks of size `N`. Each block is encrypted, and the results are concatenated.

    Args:
        msg (np.ndarray | string): Array of message data to encrypt, with values as integers.

    Returns:
        np.ndarray: Flattened array of encrypted message coefficients.

    Raises:
        Exception: If the public key parameters are not set or if encryption fails.
    """
    if isinstance(msg, str):
        msg = np.array([ord(c) % (pub_key["p"] - 1) + 1 for c in msg], dtype=np.int_)
    elif not isinstance(msg, np.ndarray):
        raise TypeError("Message must be a string or numpy array of integers.")

    print(f"Encryption funciton msg variable: {msg}")
    ntru = NTRUCipher(pub_key["N"], pub_key["p"], pub_key["q"])
    h_arr = np.array(pub_key["h"], dtype=np.int_)
    ntru.h_poly = Poly(h_arr[::-1], x).set_domain(ZZ)

    if ntru.N >= len(msg):
        msg_padded = np.pad(msg, (0, ntru.N - len(msg)), "constant")
        output = ntru.encrypt(
            Poly(msg_padded[::-1], x).set_domain(ZZ),
            random_poly(ntru.N, int(math.sqrt(ntru.q))),
        ).all_coeffs()[::-1]
    else:
        msg_arr = padding_encode(msg, ntru.N).reshape((-1, ntru.N))
        output = np.array([])
        for _, b in enumerate(msg_arr, start=1):
            next_output = ntru.encrypt(
                Poly(b[::-1], x).set_domain(ZZ),
                random_poly(ntru.N, int(math.sqrt(ntru.q))),
            ).all_coeffs()[::-1]
            if len(next_output) < ntru.N:
                next_output = np.pad(
                    next_output, (0, ntru.N - len(next_output)), "constant"
                )
            output = np.concatenate((output, next_output))

    return np.array(output).flatten()


def decrypt(msg):
    """Decrypts a message using the NTRU cryptographic scheme.

    Initializes an NTRU cipher with the private key and decrypts the input message.
    If the message length exceeds `N`, it is reshaped into blocks of size `N` and
    each block is decrypted sequentially. The resulting output is decoded to remove
    any padding that was applied during encryption.

    Args:
        msg (np.ndarray): Array of encrypted message coefficients as integers.

    Returns:
        np.ndarray: Decoded array of original message data.

    Raises:
        Exception: If the private key parameters are not set or if decryption fails.
    """
    ntru = NTRUCipher(priv_key["N"], priv_key["p"], priv_key["q"])

    f_arr = np.array(priv_key["f"], dtype=np.int_)
    f_p_arr = np.array(priv_key["f_p"], dtype=np.int_)

    ntru.f_poly = Poly(f_arr[::-1], x).set_domain(ZZ)
    ntru.f_p_poly = Poly(f_p_arr[::-1], x).set_domain(ZZ)

    if ntru.N >= len(msg):
        output = np.array(
            ntru.decrypt(Poly(msg[::-1], x).set_domain(ZZ)).all_coeffs()[::-1],
            dtype=np.int_,
        )
    else:
        msg_arr = msg.reshape((-1, ntru.N))
        output = np.array([])
        for _, b in enumerate(msg_arr, start=1):
            next_output = ntru.decrypt(Poly(b[::-1], x).set_domain(ZZ)).all_coeffs()[
                ::-1
            ]
            if len(next_output) < ntru.N:
                next_output = np.pad(
                    next_output, (0, ntru.N - len(next_output)), "constant"
                )
            output = np.concatenate((output, next_output))

        output = padding_decode(output, ntru.N)

    print(f"Decrypt Function output variable: {output}")

    if isinstance(msg, (str, bytes)) or isinstance(msg, np.ndarray):
        text = ""
        for val in output:
            if val < 0:
                val = priv_key["p"] - val
            text += chr((val - 1) % (priv_key["p"] - 1))
        return text
    return output


def main(N, p, q):
    generate_keys(N, p, q)

    print(f"Public Key: {pub_key}\nPrivate Key: {priv_key}")

    orij_text = "Hello, World!"

    enc_text = encrypt(orij_text)
    print(f"Encrypted Text: {enc_text}")
    dec_text = decrypt(enc_text)


main(701, 3, 2048)
""" Since we are using 3 as the value for p, in the encryption arrays and 
decryption arrays, the value 2 and -1 are same, because we work in mod(3) """
