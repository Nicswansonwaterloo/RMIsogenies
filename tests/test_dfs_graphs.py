from richelot_rm.genus_two_structures import GenusTwoProductStructure
from richelot_rm.product_point import ProductPoint
from richelot_rm.richelot_product_isogenies import get_arbitrary_square_example
from richelot_rm.richelot_vertex_RM import RMVertex
from tests.test_cgl_rm import (
    sqrt_two_action_on_symplectic_torsion,
)
from collections import deque
import json
from sage.all import Graph


def test_bfs_graph_traversal():
    e = 8
    p = 2**e * 5 - 1
    graph_dict = {}
    tried_initial_j_invariants = set()
    expanded = set()

    def spread_from_square():
        current_level = 0
        square = get_arbitrary_square_example(p)
        E1, E2 = square.E1, square.E2
        if E1.j_invariant() in tried_initial_j_invariants:
            raise NotImplementedError("Already tried this initial j-invariant.")

        # F_base = E1.base_field()
        # F_ext = F_base.extension(2)
        # E1, E2 = E1.change_ring(F_ext), E2.change_ring(F_ext)
        # square = GenusTwoProductStructure(E1, E2)
        working_e = e
        P2e_1, Q2e_1 = E1.torsion_basis(2 ** (working_e))
        P2e_2, Q2e_2 = E2(P2e_1), E2(Q2e_1)
        torsion_generators = [
            ProductPoint(P2e_1, E2(0)),
            ProductPoint(E1(0), P2e_2),
            ProductPoint(Q2e_1, E2(0)),
            ProductPoint(E1(0), Q2e_2),
        ]
        # initial_vertex = RMVertex(
        #     square,
        #     working_e,
        #     torsion_generators,
        #     golden_ratio_action_on_symplectic_torsion(2, working_e),
        # )
        initial_vertex = RMVertex(
            square,
            working_e,
            torsion_generators,
            sqrt_two_action_on_symplectic_torsion(2, working_e),
        )
        temp_graph_dict = {}
        queue = deque([initial_vertex])
        visited = set()
        if initial_vertex not in visited:
            visited.add(initial_vertex)

        while visited and current_level < working_e - 3:
            level_size = len(queue)
            for _ in range(level_size):
                current_vertex = queue.popleft()
                neighbors = current_vertex.get_neighbors()
                temp_graph_dict[current_vertex] = neighbors
                for neighbor in neighbors:
                    if neighbor not in visited:
                        queue.append(neighbor)
                visited.add(current_vertex)
            current_level += 1

        return temp_graph_dict, visited, initial_vertex

    for i in range(100):
        try:
            temp_graph_dict, visited, start_vert = spread_from_square()
            for k, v in temp_graph_dict.items():
                if k in graph_dict:
                    graph_dict[k] = list(set(graph_dict[k] + v))
                else:
                    graph_dict[k] = v
            expanded = expanded.union(visited)
        except (NotImplementedError, AssertionError) as ex:
            print(f"------FAIL-------")
            print(str(ex))
            continue
        print(f"Success (i = {i}): {len(expanded)}")
        tried_initial_j_invariants.add(start_vert.g2_structure[0].j_invariant())

    print(f"Number expanded: {len(expanded)}")
    G = Graph(graph_dict)
    vertex_labels = {v: v.get_type() for v in G.vertices()}
    p = G.plot(vertex_labels=vertex_labels, layout="spring")
    p.save("test_output/rm_graph/rm_full_graph.png")

    return graph_dict, start_vert


def generate_adjacency_list(graph_dict, starting_vert):
    levels = []
    adjacency_list = []
    visited = set()
    visited.add(starting_vert)
    queue = deque()
    queue.append(starting_vert)
    while queue:
        level_size = len(queue)
        current_level = []
        for _ in range(level_size):
            current_vertex = queue.popleft()
            current_level.append(current_vertex)
            neighbors = graph_dict.get(current_vertex, [])
            for neighbor in neighbors:
                adjacency_list.append((current_vertex, neighbor))
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        levels.append(current_level)

    return levels, adjacency_list


if __name__ == "__main__":
    graph_dict, starting_vert = test_bfs_graph_traversal()
    levels, adjacency_list = generate_adjacency_list(graph_dict, starting_vert)
    num_verts = 0

    def serialize_invariants(invariants):
        return hash(tuple(x.to_integer() for x in invariants))

    hashed_levels = [
        [serialize_invariants(vert.invariants) for vert in level] for level in levels
    ]
    hashed_adjacency = [
        (
            serialize_invariants(edge[0].invariants),
            serialize_invariants(edge[1].invariants),
        )
        for edge in adjacency_list
    ]
    for level in levels:
        num_verts += len(level)

    print(f"Generated graph with {num_verts} vertices and {len(adjacency_list)} edges.")
    with open("test_output/rm_graph/rm_full_graph.json", "w") as f:
        json.dump({"levels": hashed_levels, "adjacency_list": hashed_adjacency}, f)
