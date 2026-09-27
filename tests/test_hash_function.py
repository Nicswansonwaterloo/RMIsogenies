"""Digits select four non-backtracking edges; the endpoint determines the hash."""

from itertools import product

from sage.all import set_random_seed

from rm_cgl import rm_cgl, text_to_digits
from theta_rm.initialization import get_initial_vertex


set_random_seed(11)
V = get_initial_vertex(359, 3)
message = '001230'
excluded = None
for digit in message:
    edges = V.get_neighbor_edges()
    choices = [
        (W, phi, K) for W, phi, K in edges
        if (W != V if excluded is None else V.kernel_set(*K) != excluded)
    ]
    assert len(choices) == 4
    W, phi, K = choices[int(digit)]

    # The kernel of the dual isogeny is phi(A[2]).
    A = V.surface
    if A.is_product:
        basis = A.torsion_basis(2)
        points = [
            sum((a * P for a, P in zip(v, basis)), A.zero)
            for v in product(range(2), repeat=4)
        ]
    else:
        points = A.two_torsion().values()
    excluded = frozenset(phi(P) for P in points)
    assert len(excluded) == 4 and excluded == V.dual_kernel(phi)
    V = W

expected = V.rm_invariants()
for seed in (11, 22):
    set_random_seed(seed)
    assert rm_cgl(message, e=3) == expected
assert rm_cgl('', e=3) == get_initial_vertex(359, 3).rm_invariants()

assert text_to_digits('A') == '1001'
assert text_to_digits('\0A') == '00001001'
assert text_to_digits('é') == '30032221'
print('Passed: digit selection, dual kernels, seed-independent hashes and UTF-8 encoding.')
