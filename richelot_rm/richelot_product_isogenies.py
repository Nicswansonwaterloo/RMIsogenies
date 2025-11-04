from sage.all import EllipticCurve, GF, Matrix, vector, PolynomialRing
from sage.schemes.elliptic_curves.ell_finite_field import special_supersingular_curve
from sage.schemes.elliptic_curves.ell_curve_isogeny import EllipticCurveIsogeny
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

from richelot_rm.jacobian_point import JacobianPoint
from richelot_rm.product_point import ProductPoint
from richelot_rm.genus_two_structures import (
    GenusTwoProductStructure,
    GenusTwoJacobianStructure,
)


def get_square_1728_example(p):
    if p % 4 != 3:
        raise ValueError(f"p={p} must be 3 mod 4 for j=1728 to be supersingular")
    F = GF(p**2, modulus=[1, 0, 1], names="i")
    E1728 = EllipticCurve(F, [0, 0, 0, 1, 0])
    return GenusTwoProductStructure(E1728, E1728)


def get_square_0_example(p):
    if p % 3 != 2:
        raise ValueError(f"p={p} must be 2 mod 3 for j=0 to be supersingular")
    F = GF(p**2)
    E0 = EllipticCurve(F, [0, 0, 0, 0, 1]).montgomery_model()
    return GenusTwoProductStructure(E0, E0)


def get_0_and_1728_example(p):
    if p % 4 != 3:
        raise ValueError(f"p={p} must be 3 mod 4 for j=1728 to be supersingular")
    if p % 3 != 2:
        raise ValueError(f"p={p} must be 2 mod 3 for j=0 to be supersingular")
    F = GF(p**2)
    E0 = EllipticCurve(F, [0, 0, 0, 0, 1]).montgomery_model()
    E1728 = EllipticCurve(F, [0, 0, 0, 1, 0])
    return GenusTwoProductStructure(E0, E1728)


def _random_supersingular_curve(p):
    F = GF(p**2)
    E = special_supersingular_curve(F)
    K = E.random_point()
    E_rand = E.isogeny(K, algorithm="factored").codomain()

    return E_rand.montgomery_model()


def get_arbitrary_square_example(p):
    E = _random_supersingular_curve(p)
    return GenusTwoProductStructure(E, E)


def get_arbitrary_product_example(p):
    E1 = _random_supersingular_curve(p).montgomery_model()
    E2 = _random_supersingular_curve(p).montgomery_model()
    return GenusTwoProductStructure(E1, E2)


def get_0_product_example(p):
    if p % 3 != 2:
        raise ValueError(f"p={p} must be 2 mod 3 for j=0 to be supersingular")
    F = GF(p**2)
    E0 = EllipticCurve(F, [0, 0, 0, 0, 1]).montgomery_model()
    E1 = _random_supersingular_curve(p).montgomery_model()
    return GenusTwoProductStructure(E0, E1)


def get_1728_product_example(p):
    if p % 4 != 3:
        raise ValueError(f"p={p} must be 3 mod 4 for j=1728 to be supersingular")
    F = GF(p**2)
    E1728 = EllipticCurve(F, [0, 0, 0, 1, 0])
    E1 = _random_supersingular_curve(p).montgomery_model()
    return GenusTwoProductStructure(E1728, E1)


def is_2_kernel_prod(kernel):
    if not isinstance(kernel, (list, tuple)) or len(kernel) != 2:
        raise ValueError("Kernel must be a pair of CouplePoint instances.")

    gen1, gen2 = kernel
    is_isotropic = gen1.weil_pairing(gen2, 2) == 1
    is_order_2 = gen1.order() == 2 and gen2.order() == 2
    return is_isotropic and is_order_2


def is_2_kernel_diagonal(kernel):
    if not is_2_kernel_prod(kernel):
        raise ValueError("Input is not a valid 2-torsion kernel.")

    (P1, Q1), (P2, Q2) = kernel
    if P1 == 0 or P2 == 0 or Q1 == 0 or Q2 == 0:
        return True
    return False


def get_diagonal_2_isogeny(kernel):
    if not is_2_kernel_diagonal(kernel):
        raise ValueError("Input is not a diagonal 2-torsion kernel.")

    gen1, gen2 = kernel

    E1 = gen1[0].curve()
    E2 = gen1[1].curve()

    # Reformat generators to be of the form (P, 0), (0, Q)
    if gen1[0] != 0 and gen1[1] != 0:
        gen1, gen2 = gen1 + gen2, gen1
    if gen2[0] != 0 and gen2[1] != 0:
        gen1, gen2 = gen1, gen1 + gen2

    P = gen1[0] if gen1[0] != 0 else gen2[0]
    Q = gen1[1] if gen1[1] != 0 else gen2[1]

    phi1 = EllipticCurveIsogeny(E1, P)
    phi2 = EllipticCurveIsogeny(E2, Q)
    codomain = GenusTwoProductStructure(phi1.codomain(), phi2.codomain())

    def isogeny(cp_pt: ProductPoint):
        P, Q = cp_pt
        return ProductPoint(phi1(P), phi2(Q))

    return codomain, isogeny


def is_2_kernel_prod_loop(kernel):
    if not is_2_kernel_prod(kernel):
        raise ValueError("Input is not a valid 2-torsion kernel.")

    gen1, gen2 = kernel
    E1 = gen1[0].curve()
    E2 = gen1[1].curve()
    j1 = E1.j_invariant()
    j2 = E2.j_invariant()
    if j1 == j2:
        # All loops come from isomorphisms
        # if kernel is of the form (P, P), (Q, Q), then it the loop is the endomorphism:
        # [1 -1]
        # [-1 1]
        iso = E1.isomorphism_to(E2)
        if iso(gen1[0]) == gen1[1] and iso(gen2[0]) == gen2[1]:
            return True

        if j1 == 0:
            raise NotImplementedError("The case j=0 is not yet implemented.")
            # Cases: (\zeta P, P), (\zeta Q, Q) and (\zeta^2 P, P), (\zeta^2 Q, Q)
            zeta = E1.automorphisms()[2]
            P2_1 = iso(P2_1)
            Q2_1 = iso(Q2_1)
            if zeta(P2_1) == P2_2 and zeta(Q2_1) == Q2_2:
                return True
            if (zeta(zeta(P2_1)) == P2_2) and (zeta(zeta(Q2_1)) == Q2_2):
                return True
            return False

        if j1 == 1728:
            raise NotImplementedError("The case j=1728 is not yet implemented.")
            iota = E2.automorphisms()[2]  # The automorphism with iota^2 = -1
            P2_1 = iso(P2_1)
            Q2_1 = iso(Q2_1)
            return iota(P2_1) == P2_2 and iota(Q2_1) == Q2_2

    return False

def get_loop_2_isogeny(kernel):
    if not is_2_kernel_prod_loop(kernel):
        raise ValueError("Input is not an isomorphism-induced 2-torsion kernel.")
    gen1, gen2 = kernel
    E1 = gen1[0].curve()
    E2 = gen1[1].curve()
    j1 = E1.j_invariant()
    j2 = E2.j_invariant()
    if j1 == j2:
        # All loops come from isomorphisms
        # if kernel is of the form (P, P), (Q, Q), then it the loop is the endomorphism:
        # [1 -1]
        # [-1 1]
        iso = E1.isomorphism_to(E2)
        if iso(gen1[0]) == gen1[1] and iso(gen2[0]) == gen2[1]:
            def isogeny(cp_pt: ProductPoint):
                P, Q = cp_pt
                return ProductPoint(iso(P) - Q, -iso(P) + Q)
            codomain = GenusTwoProductStructure(E2, E2)
            return codomain, isogeny

    raise NotImplementedError(
        "Isomorphism-induced isogenies (LOOPS) are not yet implemented. You can check if a kernel is isomorphism-induced with is_2_kernel_isomorphism_induced(kernel) function."
    )


def is_bad_elliptic_curve_model(E):
    if E.a1() != 0 or E.a3() != 0:
        return True
    if E.a6() == 0:
        return True
    return False


def fix_curve_model(E):
    if E.a1() != 0 or E.a3() != 0:
        E_new = E.short_weierstrass_model(complete_cube=False)
        iso = E.isomorphism_to(E_new)
    else:
        E_new = E
        iso = lambda x: x

    if E_new.a6() == 0:
        Psi = WeierstrassIsomorphism(E_new, (1, -1, 0, 0))
        E_final = Psi.codomain()
        iso_final = lambda x: Psi(iso(x))
        return E_final, iso_final

    return E_new, iso


def product_to_jacobian_2_isogeny(kernel):
    if not is_2_kernel_prod(kernel):
        raise ValueError("Input is not a valid 2-torsion kernel.")

    gen1, gen2 = kernel
    E1, E2 = gen1.curves()
    if is_bad_elliptic_curve_model(E1):
        E1_iso, iso1 = fix_curve_model(E1)
    else:
        E1_iso = E1
        iso1 = lambda x: x
    if is_bad_elliptic_curve_model(E2):
        E2_iso, iso2 = fix_curve_model(E2)
    else:
        E2_iso = E2
        iso2 = lambda x: x

    Fp2 = E1_iso.base()
    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()

    P1, P2 = iso1(gen1[0]), iso2(gen1[1])
    Q1, Q2 = iso1(gen2[0]), iso2(gen2[1])

    # The roots of the cubics of E1_iso and E2_iso
    a1, a2, a3 = P1[0], Q1[0], (P1 + Q1)[0]
    b1, b2, b3 = P2[0], Q2[0], (P2 + Q2)[0]
    # Compute coefficients
    M = Matrix(Fp2, [[a1 * b1, a1, b1], [a2 * b2, a2, b2], [a3 * b3, a3, b3]])
    R, S, T = M.inverse() * vector(Fp2, [1, 1, 1])
    RD = R * M.determinant()
    da = (a1 - a2) * (a2 - a3) * (a3 - a1)
    db = (b1 - b2) * (b2 - b3) * (b3 - b1)

    s1, t1 = -da / RD, db / RD
    s2, t2 = -T / R, -S / R

    a1_t = (a1 - s2) / s1
    a2_t = (a2 - s2) / s1
    a3_t = (a3 - s2) / s1
    h = s1 * (x**2 - a1_t) * (x**2 - a2_t) * (x**2 - a3_t)
    codomain = GenusTwoJacobianStructure(h)
    J = codomain.jac

    # See https://eprint.iacr.org/2022/1283.pdf, Section 3.2.2 for derivation of formulas
    def isogeny(cp_pt: ProductPoint):
        P = iso1(cp_pt[0])
        Q = iso2(cp_pt[1])
        # The image of P
        if P != 0:
            xP, yP = P.xy()
            uP = (s1 * x**2 + s2 - xP)
            vP = Rx(yP / s1)
            div_P = J([uP, vP])
        else:
            div_P = J(0)

        # The image of Q
        if Q != 0:
            xQ, yQ = Q.xy()
            uQ = (xQ - t2) * x**2 - t1
            vQ = yQ * x**3 / t1
            div_Q = J([uQ, vQ])
        else:
            div_Q = J(0)

        return JacobianPoint(div_P + div_Q)

    return codomain, isogeny

def get_symplectic_two_torsion_prod(prod_structure: GenusTwoProductStructure):
    P1, Q1 = prod_structure.E1.torsion_basis(2)
    P2, Q2 = prod_structure.E2.torsion_basis(2)

    symplectic_basis = [
        ProductPoint(P1, prod_structure.E2(0)),
        ProductPoint(prod_structure.E1(0), P2),
        ProductPoint(Q1, prod_structure.E2(0)),
        ProductPoint(prod_structure.E1(0), Q2),
    ]

    return symplectic_basis

def compute_2_isogeny_from_product(kernel):
    if is_2_kernel_diagonal(kernel):
        return get_diagonal_2_isogeny(kernel)
    elif is_2_kernel_prod_loop(kernel):
        return get_loop_2_isogeny(kernel)
    else:
        return product_to_jacobian_2_isogeny(kernel)
