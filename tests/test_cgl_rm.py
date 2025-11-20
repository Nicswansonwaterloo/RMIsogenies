from sage.all import (
    GF,
    Matrix,
    VectorSpace,
    Integers,
    Graph,
    randint,
    set_random_seed,
    is_prime,
)
from sage.graphs.graph_latex import check_tkz_graph

from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex


def golden_ratio_action_on_symplectic_torsion(ell=2, e=1):
    """Returns the matrix of the action of the golden ratio endomorphism on the ell^e-torsion,
    with respect to the basis [ (P2_1, 0), (0, P2_2), (Q2_1, 0), (0, Q2_2) ].
    """
    Zle = Integers(ell**e)
    return Matrix(Zle, [[0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1]])


def get_cql_parameters(e):
    # Require that e \equiv 2 \pmod 4 for graph to be undirected.
    assert e % 4 == 2, "e must be congruent to 2 mod 4."
    p = 2**e * 3 - 1
    assert is_prime(p), "2**e * 3 - 1 must be prime."
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
    initial_vertex = RMVertex(
        square, e, torsion_generators, golden_ratio_action_on_symplectic_torsion(2, e)
    )
    return initial_vertex


def take_random_walk(vertex, steps, verbose=False, allow_backtrack=False):
    current_vertex = vertex
    next_vertex = None
    graph_dict = {}
    walk = []
    for step in range(steps):
        if verbose:
            print(f"Step {step}: @ vertex {current_vertex}")
        neighbors = current_vertex.get_neighbors()
        random_choice = (
            randint(0, len(neighbors) - 1)
            if allow_backtrack
            else randint(1, len(neighbors) - 1)
        )
        next_vertex = neighbors[random_choice]
        graph_dict[current_vertex] = neighbors
        current_vertex = next_vertex
        walk.append(current_vertex)

    return walk, Graph(graph_dict)


def test_random_walk(e):
    initial_vertex = get_cql_parameters(e)
    walk, G = take_random_walk(
        initial_vertex, e - 2, verbose=True, allow_backtrack=True
    )
    labels = {v: v.get_type() for v in G.vertices()}
    non_walk = [v for v in G.vertices() if v not in walk and v != initial_vertex]
    p = G.plot(
        vertex_labels=labels,
        vertex_colors={
            "#9dc3ff": walk,
            "#fadb87": non_walk,
            "#99ffa8": [initial_vertex],
        },
    )
    p.save(f"test_output/rm_graph/random_walk_e={e}.png")


def find_suitable_primes(e_start, e_end):
    suitable_primes = []
    for e in range(e_start, e_end + 1):
        p = 2**e * 3 - 1
        if is_prime(p) and e % 4 == 2:
            suitable_primes.append(e)
    return suitable_primes


def test_non_backtracking_random_walk(e):

    initial_vertex = get_cql_parameters(e)
    walk, G = take_random_walk(
        initial_vertex, e - 2, verbose=True, allow_backtrack=False
    )
    assert len(walk) == len(list(set(walk))), "Walk has backtracking!"
    labels = {v: v.get_type() for v in G.vertices()}
    non_walk = [v for v in G.vertices() if v not in walk and v != initial_vertex]
    p = G.plot(
        vertex_labels=labels,
        vertex_colors={
            "#9dc3ff": walk,
            "#fadb87": non_walk,
            "#99ffa8": [initial_vertex],
        },
    )
    p.save(f"test_output/rm_graph/random_non_backtracking_walk_e={e}.png")


if __name__ == "__main__":
    # print(find_suitable_primes(1000, 2000))
        # e = 1274
    # e = 458
    # e = 94
    # e = 34
    e = 18
    # e = 6
    check_tkz_graph()
    test_random_walk(e)
    test_non_backtracking_random_walk(e)