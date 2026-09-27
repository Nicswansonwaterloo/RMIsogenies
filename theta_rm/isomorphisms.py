"""Polarized isomorphisms between PPASs in theta coordinates."""

from sage.all import Matrix, QuadraticField, cached_function, vector

from theta_rm.theta import normalise, sqrt_fp2, theta_changes


@cached_function
def _normalizer_entries():
    # Enumerate exactly in Q(i); evaluation at a square root of -1 gives
    # the same changes of theta structure in any working characteristic.
    K = QuadraticField(-1, 'i')
    entries = tuple(tuple(M.list()) for M in theta_changes(K))
    assert len(entries) == 11520
    return entries


@cached_function
def normalizer_matrices(Fpp):
    entries = _normalizer_entries()
    i = sqrt_fp2(Fpp(-1))
    # The entries are Gaussian integers: 0, +/-1, +/-i.
    values = {c for row in entries for c in row}
    values = {c: Fpp(c[0]) + Fpp(c[1]) * i for c in values}
    return tuple(Matrix(Fpp, 4, [values[c] for c in row]) for row in entries)


@cached_function
def _normalizer_rows(Fpp):
    # The matrices share many rows. Evaluate each distinct linear form once.
    rows = {}
    indices = []
    for entries in _normalizer_entries():
        indices.append(tuple(
            rows.setdefault(entries[j:j + 4], len(rows)) for j in range(0, 16, 4)
        ))
    i = sqrt_fp2(Fpp(-1))
    rows = tuple(vector(Fpp, [Fpp(c[0]) + Fpp(c[1]) * i for c in r]) for r in rows)
    return rows, tuple(indices)


@cached_function
def theta_isomorphism_matrices(Fpp, O_A, O_B):
    """Return every normalizer matrix taking O_A to O_B projectively."""
    O_A, O_B = vector(O_A), normalise(O_B)
    rows, indices = _normalizer_rows(Fpp)
    images = [r * O_A for r in rows]
    j = next(j for j, b in enumerate(O_B) if b)
    matches = []
    for M, row_indices in zip(normalizer_matrices(Fpp), indices):
        c = images[row_indices[j]]
        if c and all(images[r] == c * b for r, b in zip(row_indices, O_B)):
            matches.append(M)
    return tuple(matches)


def theta_isomorphisms(A, B):
    """Theta-coordinate maps; f and -f agree. Use product_isomorphisms for product signs."""
    for M in theta_isomorphism_matrices(A.Fpp, A.O, B.O):
        yield lambda P, M=M: normalise(M * vector(P))
