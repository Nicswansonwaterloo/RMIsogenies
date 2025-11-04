from sage.all import GF, Matrix, VectorSpace, Integers, Graph

from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex


def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    """Returns the matrix of the action of the golden ratio endomorphism on the ell^e-torsion,
    with respect to the basis [ (P2_1, 0), (0, P2_2), (Q2_1, 0), (0, Q2_2) ].
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]])


def test_square_rm():
    p = 2**11 * 3 - 1
    e = 9
    square = get_arbitrary_square_example(p)
    E1, E2 = square.E1, square.E2
    P2e_1, Q2e_1 = E1.torsion_basis(2**e)
    P2e_2, Q2e_2 = E2(P2e_1), E2(Q2e_1)
    torsion_generators = [
        ProductPoint(P2e_1, E2(0)),
        ProductPoint(E1(0), P2e_2),
        ProductPoint(Q2e_1, E2(0)),
        ProductPoint(E1(0), Q2e_2),
    ]
    vertex = RMVertex(
        square, e, torsion_generators, golden_ratio_action_on_symplectic_torsion(2, e)
    )
    neighbors = vertex.get_neighbors()
    edges = {}
    edges[vertex] = neighbors
    G = Graph(edges)
    # Set vertex labels using get_type
    labels = {v: v.get_type() for v in G.vertices()}
    p = G.plot(vertex_labels=labels)
    p.save("square_vertex_neighbors_RM.png")

if __name__ == "__main__":
    test_square_rm()