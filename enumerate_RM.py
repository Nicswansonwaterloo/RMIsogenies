"""
This script finds all Real Multiplication (RM) endomorphisms
in the endomorphism ring of E0xE0 where E0 is a fixed supersingular elliptic
curve defined over Fp^2 with known endomorphism ring

Specifically, it constructs a lattice arising from a maximal order in a
quaternion algebra ramified at p and infinity. It then finds short vectors
in this lattice corresponding to elements with a given discriminant, trace,
and norm, and finally groups these RM candidates into equivalence classes.
"""

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


# =====================================================================
# 1. Setup
# =====================================================================
# p = 2**5 * 3**4 - 1
p = 2**216 * 3 - 1
assert is_prime(p), "p must be a prime"

D_val = 17

Fp2 = GF(p**2, "i", modulus=[1, 0, 1])
E0 = EllipticCurve(Fp2, [1, 0])
assert E0.is_supersingular(), "E0 must be supersingular"

# Define the quaternion algebra ramified at p and infinity
B = QuaternionAlgebra(QQ, -1, -p)
_, i, j, k = B.basis()

# Space of 2x2 matrices over B
H2space = MatrixSpace(B, 2, 2)

# Describe the endomorphism ring of E0 as a maximal order in B.
# Algebraic description of O0 basis
O0_basis = [B(1), i, i / 2 + j / 2, B(1 / 2) + k / 2]
O0 = B.quaternion_order(O0_basis)

# =====================================================================
# 2. Lattice Definition
# =====================================================================
lattice_basis = [
    H2space([1, 0, 0, 0]),
    H2space([0, 0, 0, 1]),
    H2space([0, O0_basis[0], O0_basis[0].conjugate(), 0]),
    H2space([0, O0_basis[1], O0_basis[1].conjugate(), 0]),
    H2space([0, O0_basis[2], O0_basis[2].conjugate(), 0]),
    H2space([0, O0_basis[3], O0_basis[3].conjugate(), 0]),
]
n = len(lattice_basis)


def herm_product(A, B_mat):
    """Compute the Hermitian inner product of two matrices."""
    return ZZ(((A * B_mat).trace()).reduced_trace())


def herm_norm(A):
    """Compute the Hermitian norm of a matrix."""
    return herm_product(A, A)


def vec_to_matrix(v):
    """Convert a vector representation back to a matrix over B."""
    return sum(QQ(coeff) * base_element for coeff, base_element in zip(v, lattice_basis))


# Construct the Gram matrix
G = matrix(QQ, n, n)
for m in range(n):
    for k in range(n):
        if m == k:
            G[m, k] = herm_product(lattice_basis[m], lattice_basis[k])

# Use LLL to find a reduced basis for the integral lattice
L = IntegralLattice(G)
L = L.LLL()


# =====================================================================
# 3. Real Multiplication Finding Algorithm
# =====================================================================
def minimal_poly_from_discriminant(D):
    """
    Given a discriminant D, returns the minimal polynomial coefficients (t, d)
    such that x^2 - t x + d = 0.
    """
    assert D % 4 in [1, 2, 3], "Discriminant must be 1, 2, or 3 mod 4"
    if D % 4 in [2, 3]:
        return 0, -D
    else:  # D % 4 == 1
        return 1, (1 - D) // 4


tr, det = minimal_poly_from_discriminant(D_val)

# For an RM element A, we want: A^2 - tr(A) + det = 0.
# The norm equation is ||A||^2 = 2 * tr^2 - 4 * det
norm_of_RM = 2 * tr**2 - 4 * det

# Find all vectors up to the required norm using LLL-reduced lattice
short_vectors = []
for batch in L.short_vectors(norm_of_RM + 1, up_to_sign_flag=True):
    # Filter non-zero length batches
    if len(batch) > 0 and any(batch[0][idx] != 0 for idx in range(n)):
        short_vectors.extend(batch)

# Filter out vectors that do not match the exact target norm
vectors_with_norm = [v for v in short_vectors if v.inner_product(v) == norm_of_RM]
print(f"Found {len(vectors_with_norm)} short vectors with the correct norm.")

# Convert vectors back to matrices and filter by trace and determinant
matrices_with_norm = [vec_to_matrix(v) for v in vectors_with_norm]
matrices_admitting_RM = []
for M in matrices_with_norm:
    if M.trace() == tr and M.determinant() == det:
        matrices_admitting_RM.append(M)

print(
    f"Found {len(matrices_admitting_RM)} matrices with the correct trace and determinant."
)


# =====================================================================
# 4. Equivalence Class Filtering
# =====================================================================
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


def conjugate_transpose(M):
    """Return the conjugate transpose of a matrix M over the quaternion algebra."""
    return M.transpose().apply_map(lambda x: x.conjugate())


def automorphisms():
    """
    Enumerate unitary actions from M_2(O_0).
    First finds norm 1 units in the maximal order O_0, then constructs
    valid diagonal and anti-diagonal matrices u.
    """
    units = []
    # Exhaustively search small coefficients for units
    for c2 in [-1, 0, 1]:
        if p * c2**2 > 1:
            continue
        for c3 in [-2, -1, 0, 1, 2]:
            if c3 != 0 and (p + 1) // 4 * c3**2 > 1:
                continue
            for c1 in [-2, -1, 0, 1, 2]:
                if c1**2 > 1:
                    continue
                for c0 in [-2, -1, 0, 1, 2]:
                    # Construct element safely within B
                    elt = B(
                        c0 * O0_basis[0]
                        + c1 * O0_basis[1]
                        + c2 * O0_basis[2]
                        + c3 * O0_basis[3]
                    )
                    if elt.reduced_norm() == 1:
                        if elt not in units:
                            units.append(elt)

    U_mats = []
    Z = B(0)
    for u1 in units:
        for u2 in units:
            # Diagonal matrices
            U_mats.append(matrix(B, 2, 2, [[u1, Z], [Z, u2]]))
            # Anti-diagonal matrices
            U_mats.append(matrix(B, 2, 2, [[Z, u1], [u2, Z]]))

    return U_mats


U_matrices = automorphisms()


def check_equivalence(A1, A2):
    """
    Check if A1 and A2 are equivalent under the action:
    A1^* u^* = u A2 for some u with u^* u = I_2.
    """
    A1_star = conjugate_transpose(A1)

    for u in U_matrices:
        u_star = conjugate_transpose(u)
        if A1_star * u_star == u * A2:
            return True, u

    return False, None


# Find equivalence class representatives (algebraic ones)
RM_equivalence_classes_reps = []
for A in matrices_admitting_RM:
    found_class = False
    for A2 in RM_equivalence_classes_reps:
        equivalent, _ = check_equivalence(A, A2)
        if equivalent:
            found_class = True
            break
    if not found_class:
        RM_equivalence_classes_reps.append(A)

# Verify that all representatives satisfy the minimal polynomial
for A in RM_equivalence_classes_reps:
    assert A**2 - tr * A + det * H2space.identity_matrix() == 0, (
        "Representative failed minimal polynomial check"
    )

print("\nRM equivalence class representatives:")
print_matrices_side_by_side(RM_equivalence_classes_reps)
