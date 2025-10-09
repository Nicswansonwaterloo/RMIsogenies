from sage.all import EllipticCurve, GF, VectorSpace, Matrix
from sage.schemes.elliptic_curves.ell_finite_field import special_supersingular_curve
from sage.schemes.elliptic_curves.ell_curve_isogeny import EllipticCurveIsogeny
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

from .dependencies.Theta_SageMath.theta_structures.couple_point import CouplePoint
from .dependencies.Castryck_Decru_SageMath.richelot_aux import FromProdToJac

# Functions to generate products of supersingular elliptic curves
# with specified or arbitrary j-invariants. The input is a F = GF(p^2) for some prime p (j = 1728 and 0 may not be supersingular over F).

def get_square_1728_example(p):
    if p % 4 != 3:
        raise ValueError(f"p={p} must be 3 mod 4 for j=1728 to be supersingular")
    F = GF(p**2, modulus=[1, 0, 1], names='i')
    E1728 = EllipticCurve(F, [0, 0, 0, 1, 0])
    return E1728, E1728

def get_square_0_example(p):
    if p % 3 != 2:
        raise ValueError(f"p={p} must be 2 mod 3 for j=0 to be supersingular")
    F = GF(p**2)
    E0 = EllipticCurve(F, [0, 0, 0, 0, 1]).montgomery_model()
    return E0, E0

def get_0_and_1728_example(p):
    if p % 4 != 3:
        raise ValueError(f"p={p} must be 3 mod 4 for j=1728 to be supersingular")
    if p % 3 != 2:
        raise ValueError(f"p={p} must be 2 mod 3 for j=0 to be supersingular")
    F = GF(p**2)
    E0 = EllipticCurve(F, [0, 0, 0, 0, 1]).montgomery_model()
    E1728 = EllipticCurve(F, [0, 0, 0, 1, 0])
    return E0, E1728

def _random_supersingular_curve(p):
    F = GF(p**2)
    E = special_supersingular_curve(F)
    K = E.random_point()
    E_rand = E.isogeny(K, algorithm="factored").codomain()

    return E_rand.montgomery_model()

def get_arbitrary_square_example(p):
    E = _random_supersingular_curve(p)
    return E, E

def get_arbitrary_product_example(p):
    E1 = _random_supersingular_curve(p).montgomery_model()
    E2 = _random_supersingular_curve(p).montgomery_model()
    return E1, E2

def get_0_product_example(p):
    if p % 3 != 2:
        raise ValueError(f"p={p} must be 2 mod 3 for j=0 to be supersingular")
    F = GF(p**2)
    E0 = EllipticCurve(F, [0, 0, 0, 0, 1]).montgomery_model()
    E1 = _random_supersingular_curve(p).montgomery_model()
    return E0, E1

def get_1728_product_example(p):
    if p % 4 != 3:
        raise ValueError(f"p={p} must be 3 mod 4 for j=1728 to be supersingular")
    F = GF(p**2)
    E1728 = EllipticCurve(F, [0, 0, 0, 1, 0])
    E1 = _random_supersingular_curve(p).montgomery_model()
    return E1728, E1

def isotropic_torsion_basis(E, N):
    """
    Returns P, Q that generate E[N]. We choose generators so that e_N(P, Q) = -1. 
    Generated naively by computing all N torsion and then checking all pairs of N-torsion points.
    This should be a deterministic function.
    """
    O = E(0)
    all_N_torsion = O.division_points(N)
    all_N_torsion = sorted(all_N_torsion)
    for i in range(len(all_N_torsion)):
        P = all_N_torsion[i]
        if P == O:
            continue
        for j in range(i+1, len(all_N_torsion)):
            Q = all_N_torsion[j]
            if Q == O:
                continue
            if P.weil_pairing(Q, N, algorithm='pari') == -1:
                return P, Q

    raise ValueError(f"Could not find basis for {N}-torsion with Weil pairing -1. Rational N torsion points: {all_N_torsion}")

def get_maximal_isotropic_subspaces(N):
    """
    Returns a list of 2-dimensional maximal isotropic subspaces of (Z/NZ)^4
    with respect to the standard symplectic form.
    """
    V = VectorSpace(GF(N), 4)
    B = Matrix(GF(N), [[0, 0, 1, 0],
                       [0, 0, 0, 1],
                       [-1, 0, 0, 0],
                       [0, -1, 0, 0]])
    subspaces = []
    for W in V.subspaces(2):
        gens = W.gens()
        u, v = gens[0], gens[1]
        if (u * B) * v == 0:
            subspaces.append(Matrix([u, v]))
    # Sort the subspaces canonically by their generators (row-wise, lex order)
    subspaces = sorted(subspaces, key=lambda mat: tuple(mat.list()))
    return subspaces

def get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, N, PN_1=None, QN_1=None, PN_2=None, QN_2=None):
    """
    Args:
        E1 (EllipticCurve): First elliptic curve over a field of characteristic != 2.
        E2 (EllipticCurve): Second elliptic curve over the same field.
        N (int): The order of the torsion subgroup (must be prime).
        PN_1, QN_1: Basis for E1[N] with Weil pairing -1 (optional).
        PN_2, QN_2: Basis for E2[N] with Weil pairing -1 (optional).

    Returns:
        list: List of tuples of CouplePoint objects generating maximal isotropic subgroups of E1[N] x E2[N].
    """
    if PN_1 is None or QN_1 is None:
        PN_1, QN_1 = isotropic_torsion_basis(E1, N)
    if PN_2 is None or QN_2 is None:
        PN_2, QN_2 = isotropic_torsion_basis(E2, N)
    symplectic_basis = [CouplePoint(PN_1, E2(0)), CouplePoint(E1(0), PN_2), CouplePoint(QN_1, E2(0)), CouplePoint(E1(0), QN_2)]
    def vec_to_point(vec):
        components = [vec[i] * symplectic_basis[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]
    subspaces = get_maximal_isotropic_subspaces(N)
    maximal_isotropic_subgroups = []
    for u, v in subspaces:
        gen1 = vec_to_point(u)
        gen2 = vec_to_point(v)
        maximal_isotropic_subgroups.append((gen1, gen2))
    return maximal_isotropic_subgroups

def is_2_kernel_diagonal(kernel):
    """
    Args:
        kernel (tuple): A tuple of two pairs, each pair containing two elements 
        (P1, P2) and (Q1, Q2) representing the generators of the kernel of a two isogeny.

    Returns:
        bool: True if the kernel generates a diagonal isogeny, False otherwise.
    """
    (P1, P2), (Q1, Q2) = kernel
    if (P2 == 0 and Q1 == 0): return True
    if (P1 == 0 and Q2 == 0): return True
    return False

def diagonal_isogeny(kernel):
    """
    Args:
        kernel (tuple): A tuple of two pairs, each pair containing two elements
                        (P1, P2) and (Q1, Q2) representing the generators of the kernel of a diagonal isogeny.
    Returns:
        tuple: A pair of two elliptic curves codomain.
                A function that takes a CouplePoint on the domain and returns its image under the isogeny.

    """
    assert is_2_kernel_diagonal(kernel), "Kernel must be diagonal."
    gen1, gen2 = kernel
    E1 = gen1[0].curve()
    E2 = gen1[1].curve()
    P1 = gen1[0] if gen1[0] != E1(0) else gen2[0]
    P2 = gen1[1] if gen1[1] != E2(0) else gen2[1]

    phi1 = EllipticCurveIsogeny(E1, P1)
    phi2 = EllipticCurveIsogeny(E2, P2)
    codomain = (phi1.codomain(), phi2.codomain())

    def isogeny(cp_pt: CouplePoint):
        P, Q = cp_pt
        return CouplePoint(phi1(P), phi2(Q))
    
    return codomain, isogeny

def is_2_kernel_isomorphism_induced(kernel):
    """
    Args:
        kernel (tuple): A tuple of two pairs of the form ((P2_1, Q2_1), (P2_2, Q2_2)) representing the generators of the kernel of a two isogeny.

    Returns:
        bool: True if the kernel generates an isomorphism-induced isogeny, meaning a loop in the graph, False otherwise.
    """
    (P2_1, P2_2), (Q2_1, Q2_2) = kernel
    E1 = P2_1.curve(); j1 = E1.j_invariant()
    E2 = P2_2.curve(); j2 = E2.j_invariant()
    if j1 == j2:
        iso = E1.isomorphism_to(E2)
        if (iso(P2_1) == P2_2 and iso(Q2_1) == Q2_2):
            return True
        
        if j1 == 0:
            raise NotImplementedError("The case j=0 is not yet implemented.")
            # Cases: (\zeta P, P), (\zeta Q, Q) and (\zeta^2 P, P), (\zeta^2 Q, Q)
            zeta = E1.automorphisms()[2]
            P2_1 = iso(P2_1); Q2_1 = iso(Q2_1)
            if (zeta(P2_1) == P2_2 and zeta(Q2_1) == Q2_2):
                return True
            if (zeta(zeta(P2_1)) == P2_2) and (zeta(zeta(Q2_1)) == Q2_2):
                return True
            return False
        
        if j1 == 1728:
            raise NotImplementedError("The case j=1728 is not yet implemented.")
            iota = E2.automorphisms()[2]  # The automorphism with iota^2 = -1
            P2_1 = iso(P2_1); Q2_1 = iso(Q2_1)
            return (iota(P2_1) == P2_2 and iota(Q2_1) == Q2_2)
        
    return False

def get_isogeny_from_product_two_kernel(kernel):
    """
    Args:
        kernel (tuple): A tuple of two pairs, each pair containing two elements
                        (P1, P2) and (Q1, Q2) representing the generators of the kernel of a two isogeny.

    Returns:
        tuple: A pair of two elliptic curves codomain or a Jacobian of a genus 2 curve.
        A function that takes a CouplePoint on the domain and returns its image under the isogeny.

    """
    if is_2_kernel_diagonal(kernel):
        return diagonal_isogeny(kernel)
    elif is_2_kernel_isomorphism_induced(kernel):
        raise NotImplementedError("Isomorphism-induced isogenies (LOOPS) are not yet implemented. You can check if a kernel is isomorphism-induced with is_2_kernel_isomorphism_induced(kernel) function.")
    else:
        (P2_1, P2_2), (Q2_1, Q2_2) = kernel
        a1, a2, a3 = P2_1[0], Q2_1[0], (P2_1 + Q2_1)[0]
        b1, b2, b3 = P2_2[0], Q2_2[0], (P2_2 + Q2_2)[0]
        E1 = P2_1.curve()
        E2 = P2_2.curve()
        Fp2 = E1.base_field()

        # The way FromProdToJac is implemented, we need to check if a certain matrix is invertible so that we can solve a system of equations. If it is not invertible, we move to a non-montgomery model of the curves and try again.

        # Compute coefficients
        M = Matrix(Fp2, [
            [a1*b1, a1, b1],
            [a2*b2, a2, b2],
            [a3*b3, a3, b3]])

        if M.determinant() != 0:
            h, _, _, _, _, isogeny = FromProdToJac(E1, E2, P2_1, Q2_1, P2_2, Q2_2, 1)
            return (h, isogeny)

        E1 = kernel[0][0].curve()
        E2 = kernel[0][1].curve()

        Psi_1 = WeierstrassIsomorphism(E1, (1, -1, 0, 0))
        E1_w = Psi_1.codomain()
        Psi_2 = WeierstrassIsomorphism(E2, (1, -1, 0, 0))
        E2_w = Psi_2.codomain()
        h, _, _, _, _, isogeny = FromProdToJac(E1_w, E2_w, Psi_1(kernel[0][0]), Psi_1(kernel[1][0]), Psi_2(kernel[0][1]), Psi_2(kernel[1][1]), 1)
        def isogeny_wrapped(cp_pt: CouplePoint):
            P, Q = cp_pt
            P = Psi_1(P)
            Q = Psi_2(Q)
            return isogeny(CouplePoint(P, Q))
        return (h, isogeny_wrapped)