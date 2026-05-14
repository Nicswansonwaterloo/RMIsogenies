from sage.all import discrete_log

from richelot_rm.product_point import ProductPoint


def get_symplectic_basis(EA, EB, M):
    B0, B1 = EA.torsion_basis(M)
    B2, B3 = EB.torsion_basis(M)
    # adjust B2 so that the pairings are the same:
    w_B0B1 = B0.weil_pairing(B1, M)
    w_B2B3 = B2.weil_pairing(B3, M)
    # we want w_B0B1 = w_B2B3, so we adjust:
    c = discrete_log(w_B0B1, w_B2B3, ord=M)
    B2 = c * B2
    
    # Standard symplectic form is
    # [ 0 0 1 0 ]
    # [ 0 0 0 1 ]
    # [ -1 0 0 0 ]
    # [ 0 -1 0 0 ]
    two_dim_basis = [
        ProductPoint(B0, EB(0)),
        ProductPoint(EA(0), B2),
        ProductPoint(B1, EB(0)),
        ProductPoint(EA(0), B3),
    ]
    
    return two_dim_basis