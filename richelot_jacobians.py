from sage.all import Matrix, GF
from dependencies.Castryck_Decru_SageMath.richelot_aux import FromJacToProd, FromJacToJac

def is_jac_kernel_split(h, kernel_generators):
    D11, D12 = kernel_generators[0]
    D21, D22 = kernel_generators[1]
    G1 = D11
    G2 = D21
    G3, _ = h.quo_rem(G1 * G2)

    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    if delta.determinant():
        # Determinant is non-zero, no splitting
        return False, None, None, None
    
    return True, G1, G2, G3

def get_isogeny_from_jacobian_two_kernel(kernel_generators):
    raise NotImplementedError("Computation of arbitrary (2,2)-isogenies from Jacobians is not yet implemented.")
    # is_split, G1, G2, G3 = is_jac_kernel_split(h, kernel_generators)
    # if is_split:
    #     isogeny, av = FromJacToProd(G1, G2, G3)
    # else:
    #     av, _, _, _, _, isogeny, _ = FromJacToJac(h, *kernel_generators[0], *kernel_generators[1], 1)
    # return av, isogeny
