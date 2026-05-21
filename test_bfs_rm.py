from collections import deque

from sage.all import Graph

from test_helpers import gen_rm_hash_prime, get_initial_vertex


def bfs_rm(initial_vertex, max_depth):
    r"""
    BFS from ``initial_vertex`` in the RM (2,2)-isogeny graph.

    INPUT:

    - ``initial_vertex`` -- starting ``RMVertex``
    - ``max_depth`` -- BFS depth cap; each step consumes one factor of 2 from
      the RM level, so vertices at depth ``d`` have ``r = initial_r - d``

    OUTPUT: pair ``(visited, adjacency)`` where ``visited`` is the set of all
    reached vertices and ``adjacency`` maps each expanded vertex to its
    neighbor list (as returned by ``get_neighbors()``)
    """
    visited = {initial_vertex}
    queue = deque([(initial_vertex, 0)])
    adjacency = {}

    while queue:
        vertex, depth = queue.popleft()
        if depth >= max_depth:
            continue

        neighbors = vertex.get_neighbors()
        adjacency[vertex] = neighbors
        print(f"Depth {depth}: {vertex} -> {len(neighbors)} neighbor(s)")

        for neighbor in neighbors:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return visited, adjacency


def test_bfs_rm():
    """Run BFS from one starting vertex in the RM graph and save a picture."""
    d = 5
    e = 4
    p, M, f = gen_rm_hash_prime(e, d)
    initial_vertex = get_initial_vertex(p, e)

    # Initial r = e + 2; stop when frontier vertices reach r = 3.
    max_depth = e

    visited, adjacency = bfs_rm(initial_vertex, max_depth)

    G = Graph(adjacency)

    labels = {v: v.get_type() for v in G.vertices()}
    non_start = [v for v in G.vertices() if v != initial_vertex]
    label = (
        rf"$K=\mathbb{{Q}}(\sqrt{{{d}}})$"
        "\n"
        rf"$p={M}\cdot 2^{{{e}}}\cdot {f} - 1$"
    )
    plot = G.plot(
        layout="spring",
        iterations=500,  # Increase iterations for a better spring layout settle
        vertex_size=200, # Adjust vertex size to fit labels better
        figsize=(10, 10), # Larger figure size distributes the nodes more widely
        vertex_labels=labels,
        vertex_colors={
            "#99ffa8": [initial_vertex],
            "#9dc3ff": non_start,
        },
        title=label,
        title_pos=(0.02, 0.98),
        fontsize=18,
    )
    filename = f"bfs_rm_d={d}_e={e}_f={f}.png"
    plot.save(filename)
    print(f"Saved {filename} with {len(visited)} vertices and {G.num_edges()} edges.")


test_bfs_rm()
