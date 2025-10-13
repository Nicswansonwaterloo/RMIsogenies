from RichelotProducts import *
from sage.all import LCM

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
    F = E1.base_field()
    P2_1, Q2_1 = isotropic_torsion_basis(E1, 2)
    assert P2_1.order() == 2
    assert Q2_1.order() == 2
    assert P2_1.weil_pairing(Q2_1, 2, algorithm='pari') == F(-1)

test_isotropic_torsion_basis()

def test_maximal_isotropic_subgroups_of_N_torsion():
    E1, E2 = get_arbitrary_product_example(2**11 * 3**5 - 1)
    maximal_isotropic_subgroups = get_maximal_isotropic_subgroups_of_N_torsion(E1, E2, 2)
    assert len(maximal_isotropic_subgroups) == 15
    seen = []
    for gen1, gen2 in maximal_isotropic_subgroups:
        # Check that the generators are valid
        assert LCM(gen1.order()[0], gen1.order()[1])  == 2
        assert LCM(gen2.order()[0], gen2.order()[1])  == 2
        assert gen1.weil_pairing(gen2, 2) == 1
        #Check there are no duplicates
        pair = (repr(gen1), repr(gen2))
        assert pair not in seen
        seen.append(pair)
    
test_maximal_isotropic_subgroups_of_N_torsion()

def test_diagonal_isogenies():
    number_of_diagonal_kernels = 0
    E1, E2 = get_arbitrary_square_example(2**11 * 3**5 - 1)
    P2_1, Q2_1 = isotropic_torsion_basis(E1, 2)
    P2_2, Q2_2 = isotropic_torsion_basis(E2, 2)
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
    
    p = 2**11 * 3**5 - 1
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