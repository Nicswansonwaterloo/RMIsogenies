"""Enumerate RM pairs (g, a) and their polarized isomorphism classes on E_0^2."""

from sage.all import (
    GF,
    EllipticCurve,
    MatrixSpace,
    QuaternionAlgebra,
    matrix,
    is_prime,
    QQ,
    ZZ,
    IntegralLattice,
)


def star(X):
    """Quaternionic conjugate transpose X^*."""
    return X.transpose().apply_map(lambda x: x.conjugate())


def trd(X):
    """Reduced trace of the matrix trace: Trd(tr(X))."""
    return X.trace().reduced_trace()


def norm_vectors(G, N):
    """All integral vectors of Gram norm N, including both signs."""
    scale = G.denominator()
    target = ZZ(scale * N)
    return IntegralLattice(scale * G).short_vectors(
        target + 1, up_to_sign_flag=False)[target]


def hermitian_basis(O0):
    """Basis B_1,...,B_6 of Sym^*(M_2(O_0))."""
    B = O0.quaternion_algebra()
    H2space = MatrixSpace(B, 2, 2)
    basis = [H2space([1, 0, 0, 0]), H2space([0, 0, 0, 1])]
    for q in O0.basis():
        basis.append(H2space([0, q, q.conjugate(), 0]))
    return basis


def get_trace_det(d):
    """Return (t, n) for K=Q(sqrt d) min polynomial x^2 - tx + n."""
    assert d > 1 and ZZ(d).is_squarefree()
    if d % 4 in (2, 3):
        return 0, -d
    return 1, (1 - d) // 4


def rm_matrices(d, O0, g):
    """Enumerate a in Sym^dagger(M_2(O_0)) with a^2 - ta + nI = 0."""
    B = O0.quaternion_algebra()
    H2space = MatrixSpace(B, 2, 2)
    I = H2space.identity_matrix()
    g_inv = g.inverse()
    basis = hermitian_basis(O0)

    # The twisted inner product on the Hermitian lattice.
    G = matrix(QQ, [[trd(g_inv * X * g_inv * Y) for Y in basis]
                    for X in basis])

    t, n = get_trace_det(d)
    N = 2 * t**2 - 4 * n

    # Sage reduces the lattice internally and returns coordinates in B_1,...,B_6.
    # Both signs are needed when t = 1.
    rms = []
    for v in norm_vectors(G, N):
        b = sum((v[k] * basis[k] for k in range(6)), H2space.zero_matrix())
        a = g_inv * b
        if a**2 - t*a + n*I == 0:
            a.set_immutable()
            rms.append(a)
    return rms


def automorphisms(O0, g):
    """Enumerate u in M_2(O_0) with u^*gu = g."""
    B = O0.quaternion_algebra()
    H2space = MatrixSpace(B, 2, 2)
    g_inv = g.inverse()
    basis = [H2space({(row, col): q})
             for row in range(2) for col in range(2) for q in O0.basis()]

    # A unitary u has Trd(tr(u^dagger u)) = 4.
    G = matrix(QQ, [[trd(g_inv * star(X) * g * Y) for Y in basis]
                    for X in basis])
    units = []
    for v in norm_vectors(G, 4):
        u = sum((v[k] * basis[k] for k in range(16)), H2space.zero_matrix())
        if star(u) * g * u == g:
            u.set_immutable()
            units.append(u)
    return units


def equivalence_classes(rms, units, g):
    """Group the RM matrices by u a_1 = a_2 u."""
    g_inv = g.inverse()
    remaining = set(rms)
    classes = []
    for a in rms:
        if a not in remaining:
            continue
        orbit = set()
        for u in units:
            b = u * a * g_inv * star(u) * g
            b.set_immutable()
            orbit.add(b)
        classes.append([b for b in rms if b in orbit])
        remaining.difference_update(orbit)
    assert not remaining
    return classes


def enumerate_rm(d, O0, g=None):
    """Return RM matrices, polarized automorphisms, and their classes."""
    B = O0.quaternion_algebra()
    g = MatrixSpace(B, 2, 2).identity_matrix() if g is None else g
    rms = rm_matrices(d, O0, g)
    units = automorphisms(O0, g)
    classes = equivalence_classes(rms, units, g)
    return rms, units, classes


def print_matrices_side_by_side(matrix_list, sep="   ", chunk_size=5):
    """Utility to print matrices neatly side-by-side."""
    if not matrix_list:
        print("[]")
        return
    for i in range(0, len(matrix_list), chunk_size):
        chunk = matrix_list[i : i + chunk_size]
        top_rows = []
        bottom_rows = []
        for M in chunk:
            strs = [[str(M[r, c]) for c in range(2)] for r in range(2)]
            c0 = max(len(strs[0][0]), len(strs[1][0]))
            c1 = max(len(strs[0][1]), len(strs[1][1]))
            top_rows.append(f"[{strs[0][0]:>{c0}}, {strs[0][1]:>{c1}}]")
            bottom_rows.append(f"[{strs[1][0]:>{c0}}, {strs[1][1]:>{c1}}]")
        print(sep.join(top_rows))
        print(sep.join(bottom_rows))
        if i + chunk_size < len(matrix_list):
            print()


if __name__ == "__main__":
    p = 31
    ds = (3, 5, 29)
    assert is_prime(p) and p % 4 == 3

    Fp2 = GF(p**2, "i", modulus=[1, 0, 1])
    E0 = EllipticCurve(Fp2, [1, 0])
    assert E0.is_supersingular()

    B = QuaternionAlgebra(QQ, -1, -p)
    _, i, j, k = B.basis()
    O0 = B.quaternion_order([B(1), i, (i + j)/2, (1 + k)/2])
    assert O0.is_maximal()
    g = MatrixSpace(B, 2, 2).identity_matrix()

    for d in ds:
        rms, units, classes = enumerate_rm(d, O0, g)
        print(f"d = {d}: {len(rms)} RM matrices, {len(units)} polarized "
              f"automorphisms, {len(classes)} "
              f"{'class' if len(classes) == 1 else 'classes'}")
        print_matrices_side_by_side([orbit[0] for orbit in classes])
        print()
