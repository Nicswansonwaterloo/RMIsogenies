# Real Multiplication CGL

This is a proof-of-concept implementation of (2,2)-isogeny graph traversal for principally polarised superspecial abelian surfaces, with tracking of a real multiplication action on the 2^r-torsion.

This code borrows snippets from
- https://github.com/GiacomoPope/Castryck-Decru-SageMath
to compute (2, 2)-isogenies using the richelot jacobian models. 

## Running

This code was tested using the SageMath version 10.9 and Python version 3.14. 

Assuming sage is in your recognised by your python interpreter, you can run examples from the repository root using `python test_cgl_rm.py` (random walk) or `python test_bfs_rm.py` (BFS traversal). These file will generate images of the vertices encountered while navigating.


Note that known errors can be thrown due to limitations in the richelot isogeny implementations. For large primes, these errors happen with negligible probability at an arbitrary vertex, but seem more likely to occur when originating from a product of elliptic curves.

## File organisation

- `enumerate_RM.py` enumerates RM endomorphisms of a starting square of elliptic curves.
- `richelot_rm/genus_two_structures.py` implements wrappers for the two surface types: `GenusTwoProductStructure` (a product of elliptic curves) and `GenusTwoJacobianStructure` (the Jacobian of a genus-2 hyperelliptic curve).
- `richelot_rm/product_point.py` implements a point on a product surface E1 × E2, with addition, scalar multiplication, order computation, and Weil pairing.
- `richelot_rm/jacobian_point.py` implements a divisor on a genus-2 Jacobian with the same interface as `ProductPoint`.
- `richelot_rm/richelot_product_isogenies.py` implements (2,2)-isogenies out of a product surface: the diagonal split, the isomorphism-induced loop, and the generic product-to-Jacobian Richelot.
- `richelot_rm/richelot_jacobian_isogeny.py` implements (2,2)-isogenies out of a Jacobian: Jacobian-to-Jacobian via the Richelot correspondence and the Jacobian-to-product split.
- `richelot_rm/richelot_vertex.py` represents a vertex in the (2,2)-isogeny graph, enumerating the 15 maximal isotropic subgroups of the 2-torsion and their corresponding neighbours, and assigns Florit–Smith type labels.
- `richelot_rm/richelot_vertex_RM.py` extends `richelot_vertex.py` to carry a 4×4 real multiplication action on the 2^r-torsion and return only RM-preserving neighbours.
- `test_helpers.py` provides shared utilities for the examples, including a symplectic 2^M-torsion basis on E1 × E2.
