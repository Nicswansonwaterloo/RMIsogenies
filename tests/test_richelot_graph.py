from sage.all import GF, Graph, PolynomialRing

from richelot_rm.genus_two_structures import GenusTwoJacobianStructure
from richelot_rm.richelot_product_isogenies import (
    get_arbitrary_product_example,
    get_arbitrary_square_example,
)
from richelot_rm.richelot_vertex import RichelotVertex


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
    assert (p ** 2 - 1) % 5 == 0
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
    assert (p**2 - 12) % 3 == 0
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



if __name__ == "__main__":
    test_square_vertex_neighbors()
    test_prod_vertex_neighbors()
    test_type_1_vertex_neighbors()
    # test_type_2_vertex_neighbors()
    test_type_3_vertex_neighbors()
    test_type_4_vertex_neighbors()
    test_type_5_vertex_neighbors()
    test_type_6_vertex_neighbors()
