from math import isqrt

from sage.all import QQ, ZZ, PolynomialRing, QuaternionAlgebra, identity_matrix, latex, matrix, vector


# Appendix: Enumerating RM. Assume p > max(d_K, 3) and p = 3 mod 4.
# Edit these radicands; any p in this range gives the same product classes.
ds = [ZZ(d) for d in (5, 29)]
discriminants = [d if d % 4 == 1 else 4*d for d in ds]
p = ZZ(max([3] + discriminants)).next_prime()
while p % 4 != 3:
    p = p.next_prime()

QA = QuaternionAlgebra(QQ, -1, -p, names=("i", "j", "k"))
i, j, k = QA.gens()
O = QA.quaternion_order([QA(1), i, (i+j)/2, (1+k)/2])
g = identity_matrix(QA, 2)


def star(a):
    return matrix(a.base_ring(), 2, 2, lambda r, s: a[s, r].conjugate())


def trd(a):
    return a[0, 0].reduced_trace() + a[1, 1].reduced_trace()


def norm_vectors(G, N):
    """All integral c with c^T G c = N; exact bounds, including both signs."""
    G = matrix(QQ, G)
    U = G.LLL_gram()  # Columns give the new basis in the original basis.
    assert abs(U.det()) == 1
    H = U.transpose() * G * U
    rank = H.nrows()

    # H = L D L^T, so q(x) = sum D[j]*(x[j]+sum L[k,j]*x[k])^2.
    L = identity_matrix(QQ, rank)
    D = []
    for j in range(rank):
        D.append(H[j, j] - sum(L[j, k]**2 * D[k] for k in range(j)))
        for r in range(j + 1, rank):
            L[r, j] = (H[r, j] - sum(L[r, k]*L[j, k]*D[k]
                                      for k in range(j))) / D[j]
    x = [ZZ(0)] * rank

    def visit(j, remaining):
        if j < 0:
            if remaining == 0:
                yield U * vector(ZZ, x)
            return
        center = -sum((L[k, j]*x[k] for k in range(j + 1, rank)), QQ(0))
        radius = isqrt(ZZ((remaining / D[j]).floor())) + 1
        for c in range(center.floor() - radius, center.ceil() + radius + 1):
            rest = remaining - D[j]*(c - center)**2
            if rest >= 0:
                x[j] = ZZ(c)
                yield from visit(j - 1, rest)

    yield from visit(rank - 1, QQ(N))


def enumerate_rm(d, O, g=None):
    """Algorithm Enumerate_RM, followed by isomorphisms of (g,a)."""
    QA = O.quaternion_algebra()
    I = identity_matrix(QA, 2)
    g = I if g is None else g
    gi = matrix(QA, [[g[1, 1], -g[0, 1]], [-g[1, 0], g[0, 0]]])
    t, n = (ZZ(1), (1-d)//4) if d % 4 == 1 else (ZZ(0), -d)
    N = 2*t*t - 4*n

    Bp = O.quaternion_algebra()
    B = [matrix(Bp, [[1, 0], [0, 0]]), matrix(Bp, [[0, 0], [0, 1]])]
    B += [matrix(Bp, [[0, q], [q.conjugate(), 0]]) for q in O.basis()]
    G = matrix(QQ, 6, 6, lambda r, s: trd(gi*B[r]*gi*B[s]))
    rms = []
    for c in norm_vectors(G, N):
        a = gi * sum((c[k]*B[k] for k in range(6)), 0*I)
        if a*a - t*a + n*I == 0:
            assert g*a == star(a)*g and all(q in O for q in a.list())
            a.set_immutable()
            rms.append(a)

    # Q_g(u) = Trd(tr(u^dagger u)) = 4; retain u^* g u = g.
    B = []
    for r in range(2):
        for s in range(2):
            for q in O.basis():
                gamma = matrix(QA, 2, 2)
                gamma[r, s] = q
                B.append(gamma)
    G = matrix(QQ, 16, 16, lambda r, s: trd(gi*star(B[r])*g*B[s]))
    units = []
    for c in norm_vectors(G, 4):
        u = sum((c[k]*B[k] for k in range(16)), 0*I)
        if star(u)*g*u == g:
            u.set_immutable()
            units.append(u)

    unseen = set(rms)
    classes = []
    for a in rms:
        if a not in unseen:
            continue
        orbit = set()
        for u in units:
            b = u*a*gi*star(u)*g  # u*a*u^{-1}, with u^{-1} = u^dagger.
            b.set_immutable()
            orbit.add(b)
        assert orbit <= unseen
        classes.append([b for b in rms if b in orbit])
        unseen.difference_update(orbit)
    assert not unseen
    return rms, units, classes


if __name__ == "__main__":
    P = PolynomialRing(QQ, "i")
    print(r"\begin{align*}")
    for index, d in enumerate(ds):
        rms, units, classes = enumerate_rm(d, O, g)
        matrices = []
        for orbit in classes:
            a = orbit[0]
            assert all(q[2] == q[3] == 0 for q in a.list())
            rows = []
            for r in range(2):
                entries = [latex(P([a[r, s][0], a[r, s][1]]))
                           .replace("i", r"\mathbf{i}") for s in range(2)]
                rows.append(" & ".join(entries))
            matrices.append(r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}")
        ending = r" \\" if index + 1 < len(ds) else ""
        print(f"d={d}: &\\quad " + r",\quad ".join(matrices) + ending)
    print(r"\end{align*}")
