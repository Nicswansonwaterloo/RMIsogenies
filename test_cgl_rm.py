from sage.all import Graph, randint

from test_helpers import gen_rm_hash_prime, get_initial_vertex


def take_random_walk(initial_vertex, e, allow_backtrack=False):
    r"""
    Take a random walk in the RM graph and return the walk and graph.

    INPUT:

    - ``initial_vertex`` -- starting ``RMVertex``
    - ``e`` -- number of 2,2-isogeny steps to take
    - ``allow_backtrack`` -- boolean (default: ``False``)

    OUTPUT: pair ``(walk, G)`` where ``walk`` is a list of vertices and ``G`` is a graph
    """
    current_vertex = initial_vertex
    # We keep track of neighbors we do not visit for testing.
    graph_dict = {}
    walk = []
    for step in range(e):
        print(f"Step {step}: @ vertex {current_vertex}")

        neighbors = current_vertex.get_neighbors()
        # neighbors are returned in a deterministic order with the backtracking neighbor first.
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


def test_non_backtracking_random_walk_cgl():
    """Draw a non-backtracking random walk picture for the fixed test case."""
    # form parameters for hash function.
    d = 5
    e = 256
    p, M, f = gen_rm_hash_prime(e, d)
    initial_vertex = get_initial_vertex(p, e)

    # take random walk.
    walk, G = take_random_walk(initial_vertex, e, allow_backtrack=False)

    # Nicely print the walk.
    labels = {v: v.get_type() for v in G.vertices()}
    non_walk = [v for v in G.vertices() if v not in walk and v != initial_vertex]
    label = (
        rf"$K=\mathbb{{Q}}(\sqrt{{{d}}})$" "\n" rf"$p={M}\cdot 2^{{{e}}}\cdot {f} - 1$"
    )
    plot = G.plot(
        layout="spring",
        iterations=500,  # Increase iterations for a better spring layout settle
        vertex_size=200,  # Adjust vertex size to fit labels better
        figsize=(10, 10),  # Larger figure size distributes the nodes more widely
        vertex_labels=labels,
        vertex_colors={
            "#9dc3ff": walk,
            "#fadb87": non_walk,
            "#99ffa8": [initial_vertex],
        },
        title=label,
        title_pos=(0.02, 0.98),
        fontsize=18,
    )
    plot.save(f"random_non_backtracking_walk_d={d}_e={e}_f={f}.png")


test_non_backtracking_random_walk_cgl()
