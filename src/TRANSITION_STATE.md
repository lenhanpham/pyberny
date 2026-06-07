# Transition-State Search (P-RFO)

`Berny(..., transition_state=True)` searches for first-order saddle points with
**partitioned rational function optimization** (P-RFO) and eigenvector
following, after Banerjee, Adams, Simons & Shepard, *J. Phys. Chem.* **89**, 52
(1985), with the Bofill/Anglada-Bofill quasi-Newton Hessian update flowchart
(*J. Comput. Chem.* **19**, 349 (1998)).

## Usage

```python
from berny import Berny, geomlib

opt = Berny(geom, transition_state=True)
for g in opt:
    energy, gradients = solver(g)
    opt.send((energy, gradients))
```

## Algorithm (`prfo_step`)

In the (projected) internal-coordinate Hessian eigenbasis each cycle:

* Diagonalize `H` into `(b_i, v_i)`; project the gradient `F_i = v_i . g`.
* The **reaction-coordinate** mode is the eigenvector of maximum overlap with
  the mode tracked on the previous cycle (lowest mode on the first cycle). It is
  shifted by the **upper** RFO root `nu_p = (b_p + sqrt(b_p^2 + 4 F_p^2)) / 2`,
  so the step climbs it.
* Every other mode shares the **lower** RFO root `nu_n`, the lowest root of the
  secular equation `nu = sum_{i != p} F_i^2 / (nu - b_i)`, so those modes are
  minimized.
* The component along mode `i` is `h_i = F_i / (nu - b_i)`; the assembled step is
  scaled uniformly into the trust radius.

The reaction-coordinate eigenvector (`BernyState.ts_mode`) is carried across
cycles, so the search keeps climbing the same physical mode.

## Hessian handling

* The model Hessian guess is seeded with one negative-curvature mode
  (`init_ts_hessian`) so P-RFO has a well-defined ascent direction.
* Updates use the flowchart `update_hessian_ts` (SR1 / BFGS / PSB selected from
  normalized cosine criteria), preserving the saddle-point signature; no linear
  search is performed (a saddle step legitimately raises the energy along the
  reaction coordinate).
* The projected Hessian is **not** forced positive-definite (the minimization
  path has no such step); a periodic spectral correction
  (`ensure_negative_eigenvalues`) restores exactly one negative eigenvalue.
* Convergence additionally requires exactly one negative Hessian eigenvalue.

## Relationship to the minimizer

The minimization path (`quadratic_step`) is unchanged. TS support is additive:
a new `transition_state` flag on `BernyParams`, a `ts_mode` field on
`BernyState`, the `prfo_step` / `update_hessian_ts` / `ensure_negative_eigenvalues`
/ `init_ts_hessian` functions, and TS branches in `Berny.send`.

## Tests

`tests/test_transition_state.py` unit-tests `prfo_step` on diagonal and rotated
quadratic saddles and drives the full optimizer to a diatomic bond maximum.
Run with `pytest` (requires numpy).
