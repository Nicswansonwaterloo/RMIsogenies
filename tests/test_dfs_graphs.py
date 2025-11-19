from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex
from tests.test_cgl_rm import golden_ratio_action_on_symplectic_torsion
from collections import deque
from sage.all import Graph

def test_dfs_graph_traversal():
    e = 7
    p = 2**e - 1
    # number of vertices:
    # (1/2880)(p − 1)(p^2 − 35p + 346) + \epsilon is about 526, enough to fit in a picture
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
    graph_dict = {}
    initial_vertex = RMVertex(
        square, e, torsion_generators, golden_ratio_action_on_symplectic_torsion(2, e)
    )
    current_level = 0
    queue = deque([initial_vertex])
    visited = {initial_vertex}
    graph_dict[initial_vertex] = []

    while queue and current_level < e - 2:
        print(f"Exploring level {current_level} with {len(queue)} vertices.")
        level_size = len(queue)
        for _ in range(level_size):
            current_vertex = queue.popleft()
            neighbors = current_vertex.get_neighbors()
            graph_dict[current_vertex] = neighbors
            for neighbor in neighbors:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
                    graph_dict[neighbor] = []
        current_level += 1

    # G = Graph(graph_dict)
    # vertex_labels = {v: v.get_type() for v in G.vertices()}
    # p = G.plot(vertex_labels=vertex_labels, layout='spring')
    # p.save("test_output/richelot_graph/rm_richelot_dfs_traversal.png")
    
if __name__ == "__main__":
    for _ in range(10):
        test_dfs_graph_traversal()