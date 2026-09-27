# RM-CGL and RMGraph traversal.

This is a proof-of-concept implementation of the RM-CGL hash function on $\mathcal{G}_K(p,\ell)$, with $K = \mathbb Q(\sqrt 5)$ and $\ell=2$.

The implementation adapts portions of the [Dartois–Maino–Pope–Robert theta-isogeny code](https://github.com/ThetaIsogenies/two-isogenies/tree/main/Theta-SageMath), under its [MIT license](vendors/two_isogenies/LICENSE).

## Running
Tested using Sage 10.10.beta10. Note: due to issues [#37109](https://github.com/sagemath/sage/issues/37109) and [#42349](https://github.com/sagemath/sage/issues/42349), it may be important to use the most recent version $>10.9$ for the corrected hyperelliptic Jacobian arithmetic.

The main two running files are `bfs_traversal.py` and `rm_cgl.py`. Given sage is in your python enviroment, these can be run from the root directory:

```sh
python bfs_traversal.py
python rm_cgl.py "RM-KLPT: Algebraic Pathfinding for Superspecial Abelian Surfaces with Maximal Real Multiplication"
```

Both scripts work over $\mathbb{F}_{p^2}$, choosing $p = 2^e \cdot 15f - 1$ with the smallest positive odd $f$ for which $p$ is prime.

`bfs_traversal.py` executes a breadth first search on $\mathcal{G}_K(p,\ell)$, and saves a plot image of the resulting graph as `RMbfs_dpth_<depth>_p=<p>.png`.

- `--depth`: maximum traversal depth (default: `4`).
- `--e`: exponent $e$ in the prime formula above (default: `3`).

`rm_cgl.py` encodes the given text as base-4 digits, follows the corresponding non-backtracking walk, and prints the RM invariants of the final vertex.

- `text`: text to hash, encoded as UTF-8 and parsed as a 4-ary string.
- `--e`: exponent $e$ in the prime formula above (default: `256`).

## A Note on Trivial Collisions

**This prototype should not be used for cryptographic applications.** As discussed in the accompanying paper, avoiding trivial collisions from loops and multi-edges can be avoided by taking $K=\mathbb{Q}(\sqrt d)$ with $d>64$, however this significantly impacts the performance of the implementation since we need to evaluate a $(d, d)$-endomorphism at each step. Thus, we stick with $d=5$ and retain such collisions at a negligible fraction of vertices, including the starting vertex.

For example, one can check that the outputs of the two messages are identical: 
```
python rm_cgl.py "o"
python rm_cgl.py "0"
```
