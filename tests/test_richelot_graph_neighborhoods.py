from sage.all import GF, Graph, PolynomialRing, is_prime, Integers, Matrix

from richelot_rm.genus_two_structures import GenusTwoJacobianStructure
from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_product_isogenies import (
    get_arbitrary_product_example,
    get_arbitrary_square_example,
    get_square_0_example,
    get_square_1728_example,
    get_0_product_example,
    get_1728_product_example,
    get_0_and_1728_example,
)
from richelot_rm.richelot_vertex import RichelotVertex
from richelot_rm.richelot_vertex_RM import RMVertex


def get_type_1_vertex(p):
    Fp2 = GF(p**2)
    s, t = Fp2.random_element(), Fp2.random_element()
    while s * t * (t**2 - 1) * (s**2 - t**2) == 0:
        s = Fp2.random_element()
        t = Fp2.random_element()

    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()
    h = (x**2 - 1) * (x**2 - s**2) * (x**2 - t**2)
    g2_structure = GenusTwoJacobianStructure(h)
    return RichelotVertex(g2_structure)


def get_type_2_vertex(p):
    assert (p**2 - 1) % 5 == 0
    Fp2 = GF(p**2)
    Fp2.zeta(5)
    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()
    h = x**5 - 1
    g2_structure = GenusTwoJacobianStructure(h)
    return RichelotVertex(g2_structure)


def get_type_3_vertex(p):
    Fp2 = GF(p**2)
    u = Fp2.random_element()
    while ((1 / u) ** 2 - 1) * (u**2 - (1 / u) ** 2) == 0:
        u = Fp2.random_element()

    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()
    h = (x**2 - 1) * (x**2 - u**2) * (x**2 - (1 / u) ** 2)
    g2_structure = GenusTwoJacobianStructure(h)
    return RichelotVertex(g2_structure)


def get_type_4_vertex(p):
    assert (p**2 - 1) % 3 == 0
    Fp2 = GF(p**2)
    zeta3 = Fp2(1).nth_root(3)
    assert zeta3**3 == 1 and zeta3 != 1
    s = 0
    t = 0
    while s * t * (t**2 - 1) * (s**2 - t**2) == 0:
        v = Fp2.random_element()
        s = ((v + 1) * (v - zeta3)) / ((v - 1) * (v + zeta3))
        t = ((v + 1) * (v - zeta3**2)) / ((v - 1) * (v + zeta3**2))

    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()
    h = (x**2 - 1) * (x**2 - s**2) * (x**2 - t**2)
    g2_structure = GenusTwoJacobianStructure(h)
    return RichelotVertex(g2_structure)


def get_type_5_vertex(p):
    assert (p**2 - 1) % 6 == 0
    Fp2 = GF(p**2)
    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()
    h = x**6 + 1
    g2_structure = GenusTwoJacobianStructure(h)
    return RichelotVertex(g2_structure)


def get_type_6_vertex(p):
    assert (p**2 - 1) % 12 == 0
    Fp2 = GF(p**2)
    Rx = PolynomialRing(Fp2, name="x")
    x = Rx.gen()
    h = x**5 + x
    g2_structure = GenusTwoJacobianStructure(h)
    return RichelotVertex(g2_structure)


def get_sage_graph_neighborhood_plot(initial_vertex: RichelotVertex):
    neighbors_with_mult = initial_vertex.get_neighbors_with_multiplicities()
    neighbors = list(neighbors_with_mult.keys())

    edges = {initial_vertex: neighbors}
    G = Graph(edges)
    vertex_labels = {v: v.get_type() for v in G.vertices()}
    for v in neighbors:
        G.set_edge_label(initial_vertex, v, neighbors_with_mult[v])

    p = G.plot(vertex_labels=vertex_labels, edge_labels=True)
    return p


# These tests are for visual inspection only, generating images in test_output/
def test_square_vertex_neighbors():
    p = 2**11 * 3 - 1
    square = get_arbitrary_square_example(p)
    initial_vertex = RichelotVertex(square)
    p = get_sage_graph_neighborhood_plot(initial_vertex)
    p.save("test_output/richelot_graph/square_vertex_neighbors.png")


def test_prod_vertex_neighbors():
    p = 2**11 * 3 - 1
    prod = get_arbitrary_product_example(p)
    initial_vertex = RichelotVertex(prod)
    p = get_sage_graph_neighborhood_plot(initial_vertex)
    p.save("test_output/richelot_graph/prod_vertex_neighbors.png")


def test_square_0_vertex_neighbors():
    p = 2**11 * 3 - 1
    square_0 = get_square_0_example(p)
    initial_vertex = RichelotVertex(square_0)
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/richelot_graph/square_0_vertex_neighbors.png")

    square_1728 = get_square_1728_example(p)
    initial_vertex = RichelotVertex(square_1728)
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/richelot_graph/square_1728_vertex_neighbors.png")

    product_0_and_1728 = get_0_and_1728_example(p)
    initial_vertex = RichelotVertex(product_0_and_1728)
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/richelot_graph/product_0_and_1728_vertex_neighbors.png")

    product_0 = get_0_product_example(p)
    initial_vertex = RichelotVertex(product_0)
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/richelot_graph/product_0_vertex_neighbors.png")

    product_1728 = get_1728_product_example(p)
    initial_vertex = RichelotVertex(product_1728)
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/richelot_graph/product_1728_vertex_neighbors.png")


def test_type_1_vertex_neighbors():
    p = 2**11 * 3 - 1
    vertex = get_type_1_vertex(p)
    p = get_sage_graph_neighborhood_plot(vertex)
    p.save("test_output/richelot_graph/type_1_vertex_neighbors.png")


def test_type_2_vertex_neighbors():
    p = 2**11 * 3 - 1
    vertex = get_type_2_vertex(p)
    p = get_sage_graph_neighborhood_plot(vertex)
    p.save("test_output/richelot_graph/type_2_vertex_neighbors.png")


def test_type_3_vertex_neighbors():
    p = 2**11 * 3 - 1
    vertex = get_type_3_vertex(p)
    p = get_sage_graph_neighborhood_plot(vertex)
    p.save("test_output/richelot_graph/type_3_vertex_neighbors.png")


def test_type_4_vertex_neighbors():
    p = 2**11 * 3 - 1
    vertex = get_type_4_vertex(p)
    p = get_sage_graph_neighborhood_plot(vertex)
    p.save("test_output/richelot_graph/type_4_vertex_neighbors.png")


def test_type_5_vertex_neighbors():
    p = 2**11 * 3 - 1
    vertex = get_type_5_vertex(p)
    p = get_sage_graph_neighborhood_plot(vertex)
    p.save("test_output/richelot_graph/type_5_vertex_neighbors.png")


def test_type_6_vertex_neighbors():
    p = 2**11 * 3 - 1
    vertex = get_type_6_vertex(p)
    p = get_sage_graph_neighborhood_plot(vertex)
    p.save("test_output/richelot_graph/type_6_vertex_neighbors.png")


def test_special_automorphic_primes():
    e = 10
    p = 2**e * 3 * 5 - 1
    product_0_and_1728 = get_0_and_1728_example(p)
    initial_vertex = RichelotVertex(product_0_and_1728)
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/richelot_graph/product_0_and_1728_vertex_neighbors.png")

    E1, E2 = product_0_and_1728
    P2e_1, Q2e_1 = E1.torsion_basis(2**e)
    P2e_2, Q2e_2 = E2.torsion_basis(2**e)
    torsion_generators = [
        ProductPoint(P2e_1, E2(0)),
        ProductPoint(E1(0), P2e_2),
        ProductPoint(Q2e_1, E2(0)),
        ProductPoint(E1(0), Q2e_2),
    ]
    Zle = Integers(2**e)
    golden_ratio_action_on_symplectic_torsion = Matrix(
        Zle, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]]
    )
    initial_vertex = RMVertex(
        product_0_and_1728,
        e,
        torsion_generators,
        golden_ratio_action_on_symplectic_torsion,
    )
    P = get_sage_graph_neighborhood_plot(initial_vertex)
    P.save("test_output/rm_graph/product_0_and_1728_vertex_neighbors.png")

if __name__ == "__main__":
    # test_square_vertex_neighbors()
    # test_prod_vertex_neighbors()
    # test_type_1_vertex_neighbors()
    # # test_type_2_vertex_neighbors()
    # test_type_3_vertex_neighbors()
    # test_type_4_vertex_neighbors()
    # test_type_5_vertex_neighbors()
    # test_type_6_vertex_neighbors()
    # test_square_0_vertex_neighbors()
    test_special_automorphic_primes()
