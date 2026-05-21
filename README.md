# RMIsogeniesRichelot

SageMath research code for walking the (2,2)-isogeny graph of principally
polarized abelian surfaces, with optional tracking of a real-multiplication
action on the 2^r-torsion. Proof of concept; not for production use.

Requires SageMath >= 10.6 (Python 3.11). Run things with `sage -python` from
the repo root.

## Files

- `richelot_rm/genus_two_structures.py` — wrappers for the two surface types:
  `GenusTwoProductStructure` (a pair of elliptic curves) and
  `GenusTwoJacobianStructure` (Jacobian of a genus-2 hyperelliptic curve).
- `richelot_rm/product_point.py` — point on `E1 x E2`, with addition, scalar
  multiplication, order, Weil pairing.
- `richelot_rm/jacobian_point.py` — divisor on a genus-2 Jacobian, same
  interface as `ProductPoint`.
- `richelot_rm/richelot_product_isogenies.py` — (2,2)-isogenies out of a
  product surface: the diagonal split case, the isomorphism-induced "loop"
  case, and the generic product-to-Jacobian Richelot.
- `richelot_rm/richelot_jacobian_isogeny.py` — (2,2)-isogenies out of a
  Jacobian: Jacobian-to-Jacobian via the Richelot correspondence and the
  Jacobian-to-product split case.
- `richelot_rm/richelot_vertex.py` — vertex in the (2,2)-graph. Enumerates the
  15 maximal isotropic subspaces of the 2-torsion and the corresponding
  neighbors. Also assigns Florian–Smith type labels.
- `richelot_rm/richelot_vertex_RM.py` — same, but carries a 4x4 RM action on
  the 2^r-torsion and only returns RM-preserving neighbors.
- `test_helpers.py` — helpers for tests, including a symplectic 2^M-torsion basis on `E1 x E2`.
- `test_cgl_rm.py` — example: non-backtracking RM-preserving random walk from
  a supersingular square, saved as a PNG.
- `test_bfs_rm.py` — placeholder for the BFS traversal from a single starting
  vertex (not yet implemented).
