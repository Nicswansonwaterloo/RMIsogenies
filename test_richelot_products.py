from richelot_products import *
from sage.all import LCM

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
        PN_1, QN_1: Basis for E1[N] (optional).
        PN_2, QN_2: Basis for E2[N] (optional).

    Returns:
        list: List of tuples of CouplePoint objects generating maximal isotropic subgroups of E1[N] x E2[N].
    """
    if PN_1 is None or QN_1 is None:
        PN_1, QN_1 = E1.torsion_basis(N)
    if PN_2 is None or QN_2 is None:
        PN_2, QN_2 = E2.torsion_basis(N)
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

def test_product_creation():
    p = 2**11 * 3**5 - 1
    E1, E2 = get_arbitrary_product_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() != E2.j_invariant()

    E1, E2 = get_square_1728_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() == E2.j_invariant() == 1728
    E1, E2 = get_square_0_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() == E2.j_invariant() == 0

    E1, E2 = get_0_and_1728_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() == 0 and E2.j_invariant() == 1728

    E1, E2 = get_0_product_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() == 0 and E2.j_invariant() != 0

    E1, E2 = get_1728_product_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() == 1728 and E2.j_invariant() != 1728

    E1, E2 = get_arbitrary_square_example(p)
    assert E1.is_supersingular() and E2.is_supersingular()
    assert E1.j_invariant() == E2.j_invariant()
    assert E1.j_invariant() != 0 and E1.j_invariant() != 1728

test_product_creation()

def test_isotropic_torsion_basis():
    E1, E2 = get_arbitrary_product_example(2**11 * 3**5 - 1)
    P2_1, Q2_1 = E1.torsion_basis(2)
    assert P2_1.order() == 2
    assert Q2_1.order() == 2
    assert P2_1.weil_pairing(Q2_1, 2).multiplicative_order() == 2

    P2e_1, Q2e_1 = E1.torsion_basis(2**11)
    assert P2e_1.order() == 2**11
    assert Q2e_1.order() == 2**11
    assert P2e_1.weil_pairing(Q2e_1, 2**11).multiplicative_order() == 2**11
    

test_isotropic_torsion_basis()

def test_maximal_isotropic_subgroups_of_N_torsion():
    E1, E2 = get_arbitrary_product_example(2**11 * 3**5 - 1)
    maximal_isotropic_subgroups = get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, 2)
    assert len(maximal_isotropic_subgroups) == 15
    seen = []
    for gen1, gen2 in maximal_isotropic_subgroups:
        # Check that the generators are valid
        assert gen1.order()  == 2
        assert gen2.order()  == 2
        assert gen1.weil_pairing(gen2, 2) == 1
        #Check there are no duplicates
        pair = (repr(gen1), repr(gen2))
        assert pair not in seen
        seen.append(pair)
    
test_maximal_isotropic_subgroups_of_N_torsion()

def test_diagonal_isogenies():
    number_of_diagonal_kernels = 0
    E1, E2 = get_arbitrary_square_example(2**11 * 3**5 - 1)
    P2_1, Q2_1 = E1.torsion_basis(2)
    P2_2, Q2_2 = E2.torsion_basis(2)
    maximal_isotropic_subgroups = get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, 2, P2_1, Q2_1, P2_2, Q2_2)
    for gen1, gen2 in maximal_isotropic_subgroups:
        if is_2_kernel_diagonal((gen1, gen2)):
            number_of_diagonal_kernels += 1
            codomain, isogeny = diagonal_isogeny((gen1, gen2))
            E1prime, E2prime = codomain
            assert isogeny(gen1) == CouplePoint(E1prime(0), E2prime(0))
            assert isogeny(gen2) == CouplePoint(E1prime(0), E2prime(0))
            assert isogeny(gen1 + gen2) == CouplePoint(E1prime(0), E2prime(0))
            P = E1.random_point()
            Q = E2.random_point()
            cp_pt = CouplePoint(P, Q)
            if 2 * cp_pt != CouplePoint(E1(0), E2(0)):
                cp_image = isogeny(cp_pt)
                assert 2 * cp_image != CouplePoint(E1prime(0), E2prime(0))

    assert number_of_diagonal_kernels == 9

test_diagonal_isogenies()

def test_isomorphism_induced_product_loops():
    #Products should never have isomorphism-induced isogenies
    p = 2**11 * 3**5 - 1
    for _ in range(3):
        E1, E2 = get_arbitrary_product_example(p)
        if E1.j_invariant() == E2.j_invariant():
            continue
        maximal_isotropic_subgroups = get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, 2)
        for gen1, gen2 in maximal_isotropic_subgroups:
            assert not is_2_kernel_isomorphism_induced((gen1, gen2))

    #Arbitrary squares have [1 1; -1 1] isogeny
    for _ in range(3):
        E1, E2 = get_arbitrary_square_example(p)
        if E1.j_invariant() == 0 or E1.j_invariant() == 1728:
            continue
        maximal_isotropic_subgroups = get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, 2)
        number_of_isomorphism_induced_kernels = 0
        for gen1, gen2 in maximal_isotropic_subgroups:
            if is_2_kernel_isomorphism_induced((gen1, gen2)):
                number_of_isomorphism_induced_kernels += 1
        assert number_of_isomorphism_induced_kernels == 1

test_isomorphism_induced_product_loops()

def test_isogeny_from_product_two_kernel():
    def count_loops_products_squares_jacobians(E1, E2):
        maximal_isotropic_subgroups = get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, 2)
        loops = 0
        products = 0
        squares = 0
        jacobians = 0
        for gen1, gen2 in maximal_isotropic_subgroups:
            if gen1[0] == E1(0) and gen2[0] == E1(0):
                print("Problematic kernel:", gen1, gen2)
            try:
                codomain, isogeny = get_isogeny_from_product_two_kernel((gen1, gen2))
            except NotImplementedError:
                loops += 1
                continue
            if isinstance(codomain, tuple):
                j1 = codomain[0].j_invariant()
                j2 = codomain[1].j_invariant()
                if j1 == j2:
                    squares += 1
                else:
                    products += 1
            else:
                jacobians += 1
        return loops, products, squares, jacobians
    
    p = 2**11 * 3 - 1
    E1, E2 = get_arbitrary_product_example(p)
    loops, products, squares, jacobians = count_loops_products_squares_jacobians(E1, E2)
    # From Florit and Smith Paper
    assert 9 == products
    assert 0 == loops
    assert 6 == jacobians
    assert 0 == squares

    E1, E2 = get_arbitrary_square_example(p)
    loops, products, squares, jacobians = count_loops_products_squares_jacobians(E1, E2)

    # From Florit and Smith Paper
    assert 6 == products
    assert 1 == loops
    assert 5 == jacobians
    assert 3 == squares

    E1, E2 = get_0_and_1728_example(p)
    loops, products, squares, jacobians = count_loops_products_squares_jacobians(E1, E2)
    # From Florit and Smith Paper
    assert 9 == products
    assert 0 == loops
    assert 6 == jacobians
    assert 0 == squares


test_isogeny_from_product_two_kernel()