from sage.all import GF, Graph

from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex import RichelotVertex

def test_square_vertex_neighbors():
    p = 2**11 * 3  - 1
    square = get_arbitrary_square_example(p)
    vertex = RichelotVertex(square)
    neighbors = vertex.get_neighbors()
    edges = {}
    edges[vertex] = neighbors
    G = Graph(edges)
    # Set vertex labels using get_type
    labels = {v: v.get_type() for v in G.vertices()}
    p = G.plot(vertex_labels=labels)
    p.save("square_vertex_neighbors.png")

if __name__ == "__main__":
    test_square_vertex_neighbors()