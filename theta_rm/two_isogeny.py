from sage.all import Matrix, vector

from vendors.two_isogenies.batched_inversion import batched_inversion

from theta_rm.theta import (
    INDICES,
    torsion_label,
    ZERO_LABEL,
    ThetaSurface,
    hadamard,
    normalise,
    sqrt_fp2,
    normalizer_generators,
)


def change_coordinates_for_kernel(A, kernel):
    """Send the supplied subgroup to {(a, 0): a in F_2^2}, whose theta action
    changes coordinate signs: x_k -> (-1)^(a.k) x_k.
    """
    two_torsion_points = {P: v for v, P in A.two_torsion().items()}
    if A.is_product:
        raise ValueError('Product kernels require elliptic curve arithmetic.')
    try:
        (a, b), (c, d) = [two_torsion_points[normalise(P)] for P in kernel]
    except (KeyError, ValueError):
        raise ValueError('Supply two rational two-torsion generators.') from None
    if (
        (a, b) == ZERO_LABEL
        or (c, d) == ZERO_LABEL
        or (a, b) == (c, d)
        or (a.dot_product(d) + b.dot_product(c))
    ):
        raise ValueError('The kernel must be a maximal isotropic subgroup of A[2].')
    K = frozenset(((a, b), (c, d), torsion_label(a + c, b + d)))
    sign_subgroup = frozenset((a, INDICES[0]) for a in INDICES[1:])
    # Induced symplectic actions, in the same order as normalizer_generators.
    actions = (
        lambda a, b: torsion_label((b[0], a[1]), (a[0], b[1])),
        lambda a, b: torsion_label((a[0], b[1]), (b[0], a[1])),
        lambda a, b: torsion_label((a[0] + b[0], a[1]), b),
        lambda a, b: torsion_label((a[0], a[1] + b[1]), b),
        lambda a, b: torsion_label((a[0], a[1] + a[0]), (b[0] + b[1], b[1])),
    )
    generators = normalizer_generators(A.Fpp)
    todo = [(K, Matrix.identity(A.Fpp, 4))]
    seen = {K}
    for L, M in todo:
        if L == sign_subgroup:
            return M
        for G, f in zip(generators, actions):
            N = frozenset(f(*v) for v in L)
            if N not in seen:
                seen.add(N)
                todo.append((N, G * M))
    raise AssertionError(
        'The changes of theta structure act transitively on the 15 maximal isotropic subgroups.'
    )


class ThetaTwoIsogeny:
    def __init__(self, A, kernel):
        self.domain = A
        self.M = change_coordinates_for_kernel(A, kernel)
        O = self.M * vector(A.O)
        H = hadamard([c * c for c in O])
        if not all(H):
            raise ValueError('A gluing quotient requires componentwise elliptic curve arithmetic.')
        AA, BB, CC, DD = H
        AA_inv, BB_inv, CC_inv, DD_inv = batched_inversion(AA, BB, CC, DD)

        B = sqrt_fp2(BB * AA_inv)
        C = sqrt_fp2(CC * AA_inv)
        D = sqrt_fp2(DD * AA_inv)
        B_inv = AA * BB_inv * B
        C_inv = AA * CC_inv * C
        D_inv = AA * DD_inv * D
        self._precomputation = (B_inv, C_inv, D_inv)

        self.roots = (A.Fpp(1), B, C, D)
        self.codomain = ThetaSurface(hadamard(self.roots))
        assert self(A.O) == self.codomain.O
        assert all(self(P) == self.codomain.O for P in kernel)

    def __call__(self, P):
        X = self.M * vector(P)
        xx, yy, zz, tt = hadamard([c * c for c in X])
        Bi, Ci, Di = self._precomputation

        yy = yy * Bi
        zz = zz * Ci
        tt = tt * Di

        return normalise(hadamard((xx, yy, zz, tt)))
