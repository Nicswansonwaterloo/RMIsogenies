from sage.all import Integers, Matrix, is_prime, kronecker

from helpers import get_symplectic_basis
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex


def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    r"""
    Return the action of the golden ratio on ell^e-torsion in a fixed basis.

    INPUT:

    - ``ell`` -- prime (default: ``2``)
    - ``e`` -- positive integer exponent (default: ``1``)

    OUTPUT: a 4x4 matrix over ``\ZZ / ell^e\ZZ``
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]])


def get_initial_vertex(p, e):
    r"""
    Return an initial RM vertex for a level $2^e$ walk and prime ``p``.

    INPUT:

    - ``p`` -- prime for the base field
    - ``e`` -- positive integer exponent

    OUTPUT: an ``RMVertex``
    """
    r = e + 2
    M = 2**r
    assert (p + 1) % M == 0, "p + 1 must be divisible by 2^(e + 2)."
    square = get_arbitrary_square_example(p)
    E1, E2 = square.E1, square.E2
    torsion_generators = get_symplectic_basis(E1, E2, M)
    rm_action = golden_ratio_action_on_symplectic_torsion(2, r)

    return RMVertex(square, r, torsion_generators, rm_action)


def gen_rm_hash_prime(e, d):
    r"""
    Generate a prime of the form $p = 4 * 2^e * f - 1$ with RM conditions.

    INPUT:

    - ``e`` -- positive integer exponent (length of walk)
    - ``d`` -- quadratic discriminant parameter for the RM

    OUTPUT: pair ``(p, M, f)`` with $p = M * 2^e * f - 1$
    """
    ell = 2
    f = 1
    # M chosen to be larges power of 2 such that m^2 > 2 + 2d.
    M = 1
    while M**2 <= 2 + 2 * d:
        M *= 2
    
    # Search for f until we get a prime p = M * 2^e * f - 1 with the right conditions
    p = M * 2**e * f - 1
    while (
        p % 4 != 3
        or not is_prime(p)
        or kronecker(d, p) != -1
        or kronecker(d, ell) != -1
    ):
        f += 1
        p = M * 2**e * f - 1
    return p, M, f
