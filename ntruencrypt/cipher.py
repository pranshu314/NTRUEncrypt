from .utils import *
from sympy.abc import x
from sympy import ZZ, Poly


class NTRUCipher:
    """NTRU Cryptographic Cipher for Polynomial-based Public Key Encryption.

    This class implements key generation, encryption, and decryption
    for the NTRU encryption scheme over a polynomial ring. The scheme
    uses specified parameters `N`, `p`, and `q` to define the polynomial
    ring, modulus values, and security properties of the system.

    Attributes:
        N (int): Degree of the polynomial ring modulus `R_poly`, typically of the form `x^N - 1`.
        p (int): Small modulus for the finite field, usually a small prime.
        q (int): Large modulus for the finite field, significantly larger than `p`.
        f_poly (Poly): Private key polynomial, which must be invertible modulo `p` and `q`.
        g_poly (Poly): Randomly generated polynomial used in the public key calculation.
        h_poly (Poly): Public key polynomial calculated from `f_poly` and `g_poly`.
        f_p_poly (Poly): Inverse of `f_poly` modulo `p`.
        f_q_poly (Poly): Inverse of `f_poly` modulo `q`.
        R_poly (Poly): Polynomial ring modulus, typically `x^N - 1`.

    Methods:
        __init__(N, p, q):
            Initializes an NTRU cipher with specified parameters and sets up the polynomial ring.
        generate_random_keys():
            Generates and sets random private and public keys for the NTRU cipher.
        generate_public_key(f_poly, g_poly):
            Computes the public key `h_poly` using the given private key `f_poly` and random polynomial `g_poly`.
        encrypt(msg_poly, rand_poly):
            Encrypts a message polynomial `msg_poly` using a random polynomial `rand_poly` and the public key.
        decrypt(msg_poly):
            Decrypts a message polynomial `msg_poly` using the private key to retrieve the original message.
    """

    N = 0
    p = 0
    q = 0
    f_poly = None
    g_poly = None
    h_poly = None
    f_p_poly = None
    f_q_poly = None
    R_poly = None

    def __init__(self, N, p, q):
        """Initialize the NTRU cipher with specified parameters.

        Sets the polynomial ring modulus and initializes the cipher parameters.

        Args:
            N (int): Degree for the polynomial ring modulus `R_poly = x^N - 1`.
            p (int): Small prime modulus.
            q (int): Larger modulus, usually much greater than `p`.
        """
        self.N = N
        self.p = p
        self.q = q
        self.R_poly = Poly(x**N - 1, x).set_domain(ZZ)

    def generate_random_keys(self):
        """Generate and set random private and public keys.

        Attempts to generate a random invertible polynomial `f_poly` and computes the public key `h_poly`
        using a random polynomial `g_poly`. If `f_poly` is not invertible after 100 tries, raises an exception.

        Raises:
            Exception: If `f_poly` cannot be generated as invertible after 20 attempts.
        """
        g_poly = random_poly(self.N, int(math.sqrt(self.q)))

        tries = 100
        while tries > 0 and (self.h_poly is None):
            f_poly = random_poly(self.N, self.N // 3, neg_ones_diff=-1)
            try:
                self.generate_public_key(f_poly, g_poly)
            except Exception as ex:
                print(f"[INFO] Try {tries} falied with Exception: {ex}")
                tries -= 1
        if self.h_poly is None:
            raise Exception("Couldn't generate invertible f")

    def generate_public_key(self, f_poly, g_poly):
        """Compute the public key polynomial `h_poly`.

        Calculates the public key `h_poly` based on the given private key `f_poly` and random polynomial `g_poly`.
        Computes modular inverses of `f_poly` modulo `p` and `q` to construct `h_poly`.

        Args:
            f_poly (Poly): Private key polynomial, which must be invertible modulo `p` and `q`.
            g_poly (Poly): Random polynomial for public key generation.
        """
        self.f_poly = f_poly
        self.g_poly = g_poly
        self.f_p_poly = invert_poly(self.f_poly, self.R_poly, self.p)
        self.f_q_poly = invert_poly(self.f_poly, self.R_poly, self.q)
        p_f_q_poly = (self.p * self.f_q_poly).trunc(self.q)
        h_before_mod = (p_f_q_poly * self.g_poly).trunc(self.q)
        self.h_poly = (h_before_mod % self.R_poly).trunc(self.q)

        self._store_coeffs_as_numpy()

    def _store_coeffs_as_numpy(self):
        """Convert and store all polynomial coefficients as numpy arrays."""
        if self.h_poly is not None:
            self.h_coeffs = np.array(self.h_poly.all_coeffs()[::-1], dtype=np.int_)
        if self.f_poly is not None:
            self.f_coeffs = np.array(self.f_poly.all_coeffs()[::-1], dtype=np.int_)
        if self.f_p_poly is not None:
            self.f_p_coeffs = np.array(self.f_p_poly.all_coeffs()[::-1], dtype=np.int_)

    def encrypt(self, msg_poly, rand_poly):
        """Encrypt a message polynomial using the public key.

        Computes the encrypted polynomial by adding the product of `rand_poly` and `h_poly` (mod `q`)
        to the message polynomial `msg_poly`.

        Args:
            msg_poly (Poly): Message polynomial to encrypt.
            rand_poly (Poly): Random polynomial to ensure encryption randomness.

        Returns:
            Poly: Encrypted polynomial representation of the message.
        """
        return (
            ((rand_poly * self.h_poly).trunc(self.q) + msg_poly) % self.R_poly
        ).trunc(self.q)

    def decrypt(self, msg_poly):
        """Decrypt a message polynomial using the private key.

        Recovers the original message polynomial by applying modular operations using `f_poly` and `f_p_poly`.
        The result is computed in two steps to ensure accuracy in modular reduction.

        Args:
            msg_poly (Poly): Encrypted polynomial message to decrypt.

        Returns:
            Poly: Decrypted polynomial, representing the original message.
        """
        a_poly = ((self.f_poly * msg_poly) % self.R_poly).trunc(self.q)
        b_poly = a_poly.trunc(self.p)
        return ((self.f_p_poly * b_poly) % self.R_poly).trunc(self.p)
