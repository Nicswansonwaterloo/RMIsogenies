from sage.all import Integers, Matrix, discrete_log, is_prime, kronecker

from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_vertex_RM import RMVertex

def get_symplectic_basis_square(E, M):
    P, Q = E.torsion_basis(M)
    two_dim_basis = [
        ProductPoint(P, E(0)),
        ProductPoint(E(0), P),
        ProductPoint(Q, E(0)),
        ProductPoint(E(0), Q),
    ]

    return two_dim_basis


def get_symplectic_basis_product(EA, EB, M):
    B0, B1 = EA.torsion_basis(M)
    B2, B3 = EB.torsion_basis(M)
    # adjust B2 so that the pairings are the same:
    w_B0B1 = B0.weil_pairing(B1, M)
    w_B2B3 = B2.weil_pairing(B3, M)
    # we want w_B0B1 = w_B2B3, so we adjust:
    c = discrete_log(w_B0B1, w_B2B3, ord=M)
    B2 = c * B2

    # Standard symplectic form is
    # [ 0 0 1 0 ]
    # [ 0 0 0 1 ]
    # [ -1 0 0 0 ]
    # [ 0 -1 0 0 ]
    two_dim_basis = [
        ProductPoint(B0, EB(0)),
        ProductPoint(EA(0), B2),
        ProductPoint(B1, EB(0)),
        ProductPoint(EA(0), B3),
    ]

    return two_dim_basis


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
    torsion_generators = get_symplectic_basis_square(square.E1, M)
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
