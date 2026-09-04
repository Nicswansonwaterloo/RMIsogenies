from theta_rm.couple_point import CouplePoint
from theta_rm.theta_point import ThetaPoint
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract


def rref_zmod(mat):
    """
    Computes the (pseudo) Reduced Row Echelon Form (RREF) of a matrix over Z/MZ.
    Assumes M = 2**e for some e, and that the matrix is full rank.
    - mat is a 2x4 matrix over Z/MZ (SageMath object)

    Returns:
    - mat: a matrix of the form:
    [1, 0, *, *]
    [0, 1, *, *]
    - basis_permutation: an array of indices tracking which original columns
                        ended up as the pivots due to column swaps.
    """
    perm = [0, 1, 2, 3]
    # Find a column with a unit (coprime to 2, i.e., an odd number) in the first row
    pivot1 = next(j for j in range(4) if mat[0, j].is_unit())
    # Swap columns to bring the pivot to index 0. Keep track of the permutation of the basis vectors.
    if pivot1 != 0:
        mat.swap_columns(0, pivot1)
        perm[0], perm[pivot1] = perm[pivot1], perm[0]

    mat.rescale_row(0, 1 / mat[0, 0])  # Inversion modulo 2^e
    mat.add_multiple_of_row(1, 0, -mat[1, 0])

    pivot2 = next(j for j in range(1, 4) if mat[1, j].is_unit())
    if pivot2 != 1:
        mat.swap_columns(1, pivot2)
        perm[1], perm[pivot2] = perm[pivot2], perm[1]

    mat.rescale_row(1, 1 / mat[1, 1])
    mat.add_multiple_of_row(0, 1, -mat[0, 1])
    return mat, perm


def reorder_basis_and_sums(basis, basis_sums, perm):
    """
    basis: list of 4 elements
    basis_sums: list of 6 elements (ordered B0+B1, B0+B2, B0+B3, B1+B2, B1+B3, B2+B3)
    perm: list of 4 indices representing the new basis order
    """
    # Hardcoded map of (min_idx, max_idx) -> index in the basis_sums array
    pair_to_idx = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}
    new_basis = [basis[i] for i in perm]
    new_basis_sums = [
        basis_sums[pair_to_idx[(min(perm[i], perm[j]), max(perm[i], perm[j]))]]
        for i, j in pair_to_idx.keys()
    ]
    return new_basis, new_basis_sums

def is_permuation(matrix):
    """Check if the given 4x4 matrix over Z/2^rZ is a (scalar) permutation matrix."""
    # A permutation matrix has exactly one nonzero entry in each row and column
    for i in range(4):
        row = matrix.row(i)
        if sum(1 for v in row if v != 0) != 1:
            return False
    for j in range(4):
        col = matrix.column(j)
        if sum(1 for v in col if v != 0) != 1:
            return False

    return True


class ThetaTorsion:
    """
    Wrapper class for the data that defines a 2^n-torsion subgroup of a Theta structure (or ProductThetaStructure).

    Designed to interface a _matrix_to_kernel method across both the ProductThetaStructure and the ThetaStructure. In the Product case, this is very straightforward.

    Due to differential addion needing P + Q, this data defining the torsion subgroup contains 10 points total. That is, if B0, B1, B2, B3 are theta point representatives of a basis of the torsion, this class stores an extended basis in the following order:
    [
        B0, B1, B2, B3,
        B0 + B1,
        B0 + B2,
        B0 + B3,
        B1 + B2,
        B1 + B3,
        B2 + B3
    ]
    """

    def __init__(self, basis: list[ThetaPoint | CouplePoint], r):
        if len(basis) == 4:
            if not all(isinstance(b, CouplePoint) for b in basis):
                raise ValueError("If given 4 basis elements must be Couple points.")
            self.basis = basis
            self.is_product = True
        elif len(basis) == 10:
            if not all(isinstance(b, ThetaPoint) for b in basis):
                raise ValueError("If given 10 basis elements must be Theta points.")
            self.basis = basis[:4]
            self.basis_sums = basis[4:]
            self.is_product = False
        else:
            raise ValueError("Basis must have length 4 or 10.")
        
        self.r = r


    def change_basis_make_torsion(self, matrix):
        """Given a 4x4 change of basis matrix over 2^r, return the corresponding new basis."""
        if (
            matrix.nrows() != 4
            or matrix.ncols() != 4
            or not all(
                isinstance(v, IntegerMod_abstract) and v.modulus() == 2**self.r
                for row in matrix
                for v in row
            )
        ):
            raise ValueError(
                f"Matrix must be of size 4x4 with entries integers mod 2^r = {2**self.r}."
            )
        if is_permuation(matrix):
            print("is permuation")

        return None




    def matrix_to_kernel(self, matrix):
        """Given a 2x4 matrix over 2^r, return the corresponding 2 kernel points in a tuple by interpretting column vecs as points."""
        if (
            matrix.nrows() != 4
            or matrix.ncols() != 2
            or not all(
                isinstance(v, IntegerMod_abstract) and v.modulus() == 2**self.r
                for row in matrix
                for v in row
            )
        ):
            print(matrix)
            print(matrix.parent())
            raise ValueError(
                f"Matrix must be of size 2x4 with entries integers mod 2^r = {2**self.r}."
            )

        if self.is_product:
            col1 = matrix.column(0)
            col2 = matrix.column(1)
            components1 = [col1[i].lift_centered() * self.basis[i] for i in range(4)]
            components2 = [col2[i].lift_centered() * self.basis[i] for i in range(4)]
            K1 = components1[0] + components1[1] + components1[2] + components1[3]
            K2 = components2[0] + components2[1] + components2[2] + components2[3]
            return K1, K2

        reduced_matrix, basis_perm = self.rref_zmod(matrix.transpose())
        print(reduced_matrix)
        basis, basis_sums = reorder_basis_and_sums(
            self.basis, self.basis_sums, basis_perm
        )
        v02, v03 = int(reduced_matrix[0, 2]), int(reduced_matrix[0, 3])
        v12, v13 = int(reduced_matrix[1, 2]), int(reduced_matrix[1, 3])

        # v02B2 = tc_0.Mult(basis[2], v02)
        v02B2 = basis[2] * v02
        # v03B3 = tc_0.Mult(basis[3], v03)
        v03B3 = basis[3] * v03
        # B0_p_v02B2 = tc_0.Kxpy_xpy(v02, basis[2], basis[0], basis_sums[1])  # B0 + v02*B2
        B0_p_v02B2 = basis[0].diff_add_ladder(basis[2], basis_sums[1], v02)
        # B0_p_v03B3 = tc_0.Kxpy_xpy(v03, basis[3], basis[0], basis_sums[2])  # B0 + v03*B3
        B0_p_v03B3 = basis[0].diff_add_ladder(basis[3], basis_sums[2], v03)
        # B2_p_v03B3 = tc_0.Kxpy_xpy(v03, basis[3], basis[2], basis_sums[5])  # B2 + v03*B3
        B2_p_v03B3 = basis[2].diff_add_ladder(basis[3], basis_sums[5], v03)
        # v02B2_p_v03B3 = tc_0.Kxpy_xpy(v02, basis[2], v03B3, B2_p_v03B3)  # v02*B2 + v03*B3
        v02B2_p_v03B3 = v03B3.diff_add_ladder(
            basis[2], B2_p_v03B3, v02
        )  # v02*B2 + v03*B3
        # K1 = tc_0.Extended_Addition(
        #     B0, v02B2, v03B3, B0_p_v02B2, v02B2_p_v03B3, B0_p_v03B3
        # )  # B0 + v02*B2 + v03*B3
        K1 = basis[0].extended_addition(
            v02B2, v03B3, B0_p_v02B2, v02B2_p_v03B3, B0_p_v03B3
        )

        # v12B2 = tc_0.Mult(basis[2], v12)
        v12B2 = basis[2] * v12
        # v13B3 = tc_0.Mult(basis[3], v13)
        v13B3 = basis[3] * v13
        # B1_p_v12B2 = tc_0.Kxpy_xpy(v12, basis[2], basis[1], basis_sums[3])  # B1 + v12*B2
        B1_p_v12B2 = basis[1].diff_add_ladder(
            basis[2], basis_sums[3], v12
        )  # B1 + v12*B2
        # B1_p_v13B3 = tc_0.Kxpy_xpy(v13, basis[3], basis[1], basis_sums[4])  # B1 + v13*B3
        B1_p_v13B3 = basis[1].diff_add_ladder(
            basis[3], basis_sums[4], v13
        )  # B1 + v13*B3
        # B2_p_v13B3 = tc_0.Kxpy_xpy(v13, basis[3], basis[2], basis_sums[5])  # B2 + v13*B3
        B2_p_v13B3 = basis[2].diff_add_ladder(
            basis[3], basis_sums[5], v13
        )  # B2 + v13*B3
        # v12B2_p_v13B3 = tc_0.Kxpy_xpy(v12, basis[2], v13B3, B2_p_v13B3)  # v12*B2 + v13*B3
        v12B2_p_v13B3 = v13B3.diff_add_ladder(
            basis[2], B2_p_v13B3, v12
        )  # v12*B2 + v13*B3
        # K2 = tc_0.Extended_Addition(
        #     B1, v12B2, v13B3, B1_p_v12B2, v12B2_p_v13B3, B1_p_v13B3
        # )  # B1 + v12*B2 + v13*B3
        K2 = basis[1].extended_addition(
            v12B2, v13B3, B1_p_v12B2, v12B2_p_v13B3, B1_p_v13B3
        )
        return K1, K2

    def __getitem__(self, i):
        return self.extended_basis[i]

    def __len__(self):
        return len(self.extended_basis)

    def get_projection(self, m):
        """Returns 2^(n - m) times the basis elements, yeilding a new ThetaTorsion Object for the 2^m torsion."""
        if m < 0 or m >= self.r:
            raise ValueError("m is out of bounds.")
        new_basis = [b * (2 ** (self.r - m)) for b in self.basis]
        return ThetaTorsion(new_basis, m)


# class ThetaKernel:
#     """
#     Wrapper class for the data that defines a kernel of a  subgroup of a Theta structure.

#     Due to differential addion needing P + Q, this thetaTorsion subgroup contains 10 points total. That is, if B0, B1, B2, B3 are theta point representatives of a basis of the two-torsion, this class stores an extended basis in the following order:
#     [
#     B0, B1, B2, B3,
#     B0 + B1, B0 + B2, B0 + B3, B1 + B2, B1 + B3, B2 + B3
#     ]
#     """

#     def __init__(self, n, extended_basis: list[ThetaPoint]):
#         self.n = n
#         self.extended_basis = extended_basis

#     def __getitem__(self, i):
#         return self.extended_basis[i]

#     def __len__(self):
#         return len(self.extended_basis)
