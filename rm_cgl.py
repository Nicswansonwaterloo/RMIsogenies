from argparse import ArgumentParser
from time import perf_counter

from theta_rm.initialization import gen_rm_hash_prime, get_initial_vertex


def text_to_digits(text):
    """Encode UTF-8 bytes as four base-4 digits each, most significant first."""
    digits = []
    for byte in text.encode('utf-8'):
        for place in (64, 16, 4, 1):
            digit = (byte // place) % 4
            digits.append(str(digit))
    return ''.join(digits)


def rm_cgl(digits, e=256):
    r"""Hash quaternary digits by a non-backtracking walk for K = Q(sqrt(5))."""
    if not isinstance(digits, str) or any(d not in '0123' for d in digits):
        raise ValueError('Expected a string of digits 0 through 3.')
    if e < 3:
        raise ValueError('The prime exponent e must be at least 3.')
    p, _ = gen_rm_hash_prime(e)
    vertex = get_initial_vertex(p, e)
    dual_kernel = None
    for step, digit in enumerate(digits, start=1):
        print(f'Step {step}/{len(digits)}', flush=True)
        edges = vertex.get_neighbor_edges(dual_kernel)
        # Exclude the loop at the initial vertex.
        choices = (
            [edge for edge in edges if edge[0] != vertex]
            if dual_kernel is None else edges[1:]
        )
        assert len(choices) == 4
        next_vertex, phi, _ = choices[int(digit)]
        dual_kernel = vertex.dual_kernel(phi)
        vertex = next_vertex
    return vertex.rm_invariants()


def main():
    parser = ArgumentParser(description=main.__doc__)
    parser.add_argument(
        'text', nargs='?',
        default='RM-KLPT: Algebraic Pathfinding for Superspecial Abelian Surfaces with Maximal Real Multiplication',
        help='text to hash (encoded as UTF-8; default: %(default)s)',
    )
    parser.add_argument('--e', type=int, default=256, help='prime exponent (default: 256)')
    args = parser.parse_args()
    digits = text_to_digits(args.text)
    if args.e < 3:
        parser.error('e must be at least 3')
    start_time = perf_counter()
    kind, theta_null, kernel, action = rm_cgl(digits, args.e)
    end_time = perf_counter()
    print(f'RM canonical output ({kind}):')
    print('A[3] determining theta coords:' if kind == 'J' else 'A[3] determining product coords:')
    for i, entry in enumerate(action, start=1):
        label = (
            f'theta(P_{i}), theta(tau(P_{i})), theta(P_{i} + tau(P_{i}))'
            if kind == 'J' else f'P_{i}, tau(P_{i})'
        )
        print(f'  ({label}) = {entry}')
    print('tau kernel:')
    for i, P in enumerate(kernel, start=1):
        label = f'theta(Q_{i})' if kind == 'J' else f'Q_{i}'
        print(f'  {label} = {P}')
    print(f'theta null: {theta_null}')
    print(f'Hash completed in {end_time - start_time:.1f} seconds.')

if __name__ == '__main__':
    from sage.all import proof, set_random_seed

    set_random_seed(0)
    proof.all(False)
    main()
