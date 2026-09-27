from argparse import ArgumentParser
from collections import deque
from time import perf_counter

from sage.all import DiGraph

from theta_rm.initialization import gen_rm_hash_prime, get_initial_vertex


def bfs_rm(initial_vertex, max_depth):
    r"""Perform a breadth-first search on the RM-graph with K = Q(\sqrt{5}) and \ell = 2."""
    visited = {initial_vertex}
    queue = deque([(initial_vertex, 0)])
    adjacency = {}

    while queue:
        vertex, depth = queue.popleft()
        if depth >= max_depth:
            continue

        neighbors = [
            V
            for V, m in vertex.get_neighbors().items()
            for _ in range(m)
        ]
        adjacency[vertex] = neighbors
        neighbor_types = [V.get_type() for V in dict.fromkeys(neighbors)]
        print(f"Depth {depth}: {vertex.get_type()} -> {neighbor_types}")

        for neighbor in neighbors:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return visited, adjacency


def main():
    """Run BFS from one starting vertex in the RM graph and save a picture."""
    parser = ArgumentParser(description=main.__doc__)
    parser.add_argument(
        "--e", type=int, default=3, help="prime exponent e >= 3 (default: 3)"
    )
    parser.add_argument(
        "--depth", type=int, default=4, help="BFS depth cap (default: 4)"
    )
    args = parser.parse_args()
    if args.e < 3 or args.depth < 0:
        parser.error("e must be at least 3 and depth must be nonnegative")
    d = 5
    e = args.e
    p, f = gen_rm_hash_prime(e, d)
    initial_vertex = get_initial_vertex(p, e)

    max_depth = args.depth

    start = perf_counter()
    visited, adjacency = bfs_rm(initial_vertex, max_depth)
    print(f"Graph traversal completed in {perf_counter() - start:.1f} seconds!")

    G = DiGraph(adjacency, format="dict_of_lists", loops=True, multiedges=True)

    # relabel graph: use integer labels to avoid RM isomorphism checks while plotting.
    vertices = list(G)
    indices = {v: i for i, v in enumerate(vertices)}
    labels = {i: v.get_type() for i, v in enumerate(vertices)}
    start_index = indices.get(initial_vertex)
    G.relabel(indices)
    non_start = [i for i in G if i != start_index]

    label = (
        rf"$K=\mathbb{{Q}}(\sqrt{{{d}}})$"
        "\n"
        rf"$p=2^{{{e}}}\cdot 5\cdot 3\cdot {f} - 1$"
    )
    filename = f"RMbfs_dpth_{max_depth}_p={p}.png"
    print(f"Saving Graph image to {filename}")
    plot = G.plot(
        layout="spring",
        vertex_size=30 * 2 ** (-max(0, max_depth - 4)),
        arrowsize=0.6,
        edge_thickness=0.5,
        dist=0.015,
        figsize=(16, 16),
        vertex_labels=labels,
        label_fontsize=7,
        vertex_colors={
            "#99ffa8": [start_index] if start_index is not None else [],
            "#9dc3ff": non_start,
        },
        title=label,
        title_pos=(0.02, 0.98),
        fontsize=28,
    )
    plot.save(filename)

if __name__ == "__main__":
    main()
