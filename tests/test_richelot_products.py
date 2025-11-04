from itertools import chain, product
from richelot_rm.richelot_product_isogenies import *
from sage.all import GF, Matrix, VectorSpace, ZZ, randint


def generate_point_order_N(E, N):
    p = E.base().characteristic()
    n = (p + 1) // N
    for _ in range(1000):
        P = E.random_point()
        Q = n * P
        if Q.order() == N:
            return Q

    raise ValueError(f"Never found a point P of order N.")

def get_all_2_kernels(E1, E2, randomize_generators=False):
    symplectic_basis = get_symplectic_two_torsion_prod(GenusTwoProductStructure(E1, E2))

    V = VectorSpace(GF(2), 4)
    B = Matrix(GF(2), [[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
    subspaces = []
    for W in V.subspaces(2):
        gens = W.gens()
        u, v = gens[0], gens[1]
        if (u * B) * v == 0:
            subspaces.append(Matrix([u, v]))

    def vec_to_point(vec):
        components = [ZZ(vec[i]) * symplectic_basis[i] for i in range(4)]
        return components[0] + components[1] + components[2] + components[3]

    maximal_isotropic_subgroups = []
    for u, v in subspaces:
        gen1 = vec_to_point(u)
        gen2 = vec_to_point(v)
        if randomize_generators:
            if bool(randint(0, 1)):
                if bool(randint(0, 1)):
                    gen1 = gen1 + gen2
                else:
                    gen2 = gen1 + gen2
        kernel = ProductPoint(gen1, gen2)
        maximal_isotropic_subgroups.append(kernel)

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


def test_maximal_isotropic_subgroups_of_N_torsion():
    E1, E2 = get_arbitrary_product_example(2**11 * 3**5 - 1)
    maximal_isotropic_subgroups = get_all_2_kernels(E1, E2, randomize_generators=True)
    assert len(maximal_isotropic_subgroups) == 15
    seen = []
    for gen1, gen2 in maximal_isotropic_subgroups:
        # Check that the generators are valid
        assert gen1.order() == 2
        assert gen2.order() == 2
        assert gen1.weil_pairing(gen2, 2) == 1
        # Check there are no duplicates
        pair = (repr(gen1), repr(gen2))
        assert pair not in seen
        seen.append(pair)


def test_diagonal_isogenies():
    number_of_diagonal_kernels = 0
    E1, E2 = get_arbitrary_square_example(2**11 * 3**5 - 1)
    maximal_isotropic_subgroups = get_all_2_kernels(E1, E2, randomize_generators=True)
    for gen1, gen2 in maximal_isotropic_subgroups:
        if is_2_kernel_diagonal((gen1, gen2)):
            number_of_diagonal_kernels += 1
            codomain, isogeny = get_diagonal_2_isogeny((gen1, gen2))
            E1prime, E2prime = codomain
            assert isogeny(gen1) == ProductPoint(E1prime(0), E2prime(0))
            assert isogeny(gen2) == ProductPoint(E1prime(0), E2prime(0))
            assert isogeny(gen1 + gen2) == ProductPoint(E1prime(0), E2prime(0))
            P = E1.random_point()
            Q = E2.random_point()
            cp_pt = ProductPoint(P, Q)
            if 2 * cp_pt != ProductPoint(E1(0), E2(0)):
                cp_image = isogeny(cp_pt)
                assert 2 * cp_image != ProductPoint(E1prime(0), E2prime(0))

    assert number_of_diagonal_kernels == 9


def test_isomorphism_induced_product_loops():
    # Products should never have isomorphism-induced isogenies
    p = 2**11 * 3**5 - 1
    for _ in range(3):
        E1, E2 = get_arbitrary_product_example(p)
        if E1.j_invariant() == E2.j_invariant():
            continue
        maximal_isotropic_subgroups = get_all_2_kernels(
            E1, E2, randomize_generators=True
        )
        for gen1, gen2 in maximal_isotropic_subgroups:
            assert not is_2_kernel_prod_loop((gen1, gen2))

    # Arbitrary squares have [1 1; -1 1] isogeny
    for _ in range(3):
        E1, E2 = get_arbitrary_square_example(p)
        if E1.j_invariant() == 0 or E1.j_invariant() == 1728:
            continue
        maximal_isotropic_subgroups = get_all_2_kernels(
            E1, E2, randomize_generators=True
        )
        number_of_isomorphism_induced_kernels = 0
        for gen1, gen2 in maximal_isotropic_subgroups:
            if is_2_kernel_prod_loop((gen1, gen2)):
                number_of_isomorphism_induced_kernels += 1
        assert number_of_isomorphism_induced_kernels == 1


def test_isogeny_from_product_two_kernel():
    def count_loops_products_squares_jacobians(E1, E2):
        maximal_isotropic_subgroups = get_all_2_kernels(
            E1, E2, randomize_generators=True
        )
        loops = 0
        products = 0
        squares = 0
        jacobians = 0
        for gen1, gen2 in maximal_isotropic_subgroups:
            try:
                codomain, isogeny = compute_2_isogeny_from_product((gen1, gen2))
            except NotImplementedError as e:
                loops += 1
                continue
            if codomain.is_product:
                j1 = codomain[0].j_invariant()
                j2 = codomain[1].j_invariant()
                if j1 == j2:
                    squares += 1
                else:
                    products += 1
            else:
                # print(f"{E1.j_invariant()}, {E2.j_invariant()} -> {codomain}")
                jacobians += 1
        return loops, products, squares, jacobians

    p = 2**11 * 3 - 1
    E1, E2 = get_arbitrary_product_example(p)
    loops, products, squares, jacobians = count_loops_products_squares_jacobians(E1, E2)
    # From Florit and Smith Paper
    assert 15 == products + squares + jacobians + loops
    assert 9 == products + squares, f"Got {products} and {squares} squares"
    assert 0 == loops
    assert 6 == jacobians

    E1, E2 = get_arbitrary_square_example(p)
    loops, products, squares, jacobians = count_loops_products_squares_jacobians(E1, E2)

    # From Florit and Smith Paper
    assert 6 == products, f"Got {products} and {squares} squares"
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

def push_2e_torsion_through_M_2_isogenies(p, e, M):
    product = get_arbitrary_product_example(p)
    chain_of_products = [product]
    chain_of_isogenies = []
    for _ in range(M):
        maximal_isotropic_subgroups = get_all_2_kernels(product[0], product[1], randomize_generators=True)
        product_neighbors = []
        for gen1, gen2 in maximal_isotropic_subgroups:
            try:
                codomain, isogeny = compute_2_isogeny_from_product((gen1, gen2))
            except NotImplementedError: #loop
                continue

            if codomain.is_product:
                product_neighbors.append((codomain, isogeny))
        
        if len(product_neighbors) == 0:
            break

        random_neighbor = product_neighbors[randint(0, len(product_neighbors)-1)]
        chain_of_products.append(random_neighbor[0])
        chain_of_isogenies.append(random_neighbor[1])
        product = random_neighbor[0]

    # Pick a random 2^e torsion point on initial product
    E1, E2 = chain_of_products[0]
    T1 = generate_point_order_N(E1, 2**e)
    T2 = generate_point_order_N(E2, 2**e)
    T = ProductPoint(T1, T2)
    for isogeny in chain_of_isogenies:
        T = isogeny(T)
    return T

def test_product_chain():
    e = 11
    p = 2**e * 3 - 1
    # T has about a 1/9 chance of being in the kernel at each step So, the expected value after 9*e steps is about 0
    M = 9 * e
    images_of_T_torsions_order = []
    num_trails = 10
    for _ in range(num_trails):
        T = push_2e_torsion_through_M_2_isogenies(p, e, M)
        images_of_T_torsions_order.append(T.order())

    print(sorted(images_of_T_torsions_order))
    assert sum(images_of_T_torsions_order) / num_trails < 1 + 1e-5

if __name__ == "__main__":
    test_product_creation()
    test_maximal_isotropic_subgroups_of_N_torsion()
    test_diagonal_isogenies()
    test_isomorphism_induced_product_loops()
    for i in range(1):
        test_isogeny_from_product_two_kernel()
    test_product_chain()


