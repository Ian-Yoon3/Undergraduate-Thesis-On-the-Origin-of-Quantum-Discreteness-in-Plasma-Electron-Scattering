# p-adic Structure of the Amplituhedron

Source code for the undergraduate thesis
**On the Origin of Quantum Discreteness in Plasma Electron Scattering. A Non-Archimedean Exploration of the Amplituhedron and Its p-adic Structure**
(Sung-hyuk Yoon, Dongguk University, 2026. Supervised by Dr. Jeon ByungHoon.)

## Overview

The code computes the canonical form of the amplituhedron in two parallel tracks.
The QQ track works in exact rational arithmetic.
The QP track adds a p-adically scaled perturbation to the same input.
Both tracks are exact, so the comparison is not a comparison of floating-point error.
The difference between the two tracks is read through its p-adic valuation instead of its magnitude.

## Scope

The results do not show that a p-adic framework computes the amplitude more accurately.
They do not show that probabilistic uncertainty is reduced either.
What the code demonstrates is narrower.
The choice of number system leaves a definite and reproducible imprint on an exactly computed amplitude.

## Algorithm

The implementation follows the five-step algorithm of Section 3.3.3 of the thesis.

1. Generate the moment-curve matrix Z.
2. Test the positivity condition on the ordered maximal minors of Z.
3. Decompose the amplituhedron into BCFW cells.
4. Evaluate the canonical form of each cell.
5. Sum the cell contributions into the total amplitude.

## Why plain Python instead of SageMath

The computation only needs exact rational arithmetic and p-adic valuations.
Python's standard `fractions` module provides the first.
The second is implemented directly in a few lines.
The SageMath objects `QQ`, `Qp` and `matrix(QQ, ...)` are therefore replaced by small standalone classes.
The code runs with NumPy and Matplotlib only.

## Files

| File | Content |
|---|---|
| `experiment1_standalone.py` | Core implementation and Experiment 1 (QQ vs QP as the particle number n grows) |
| `experiment2_matrix.py` | Experiment 2 (grid over precision and n_step). Imports the functions of Experiment 1 |
| `requirements.txt` | Package versions the code was tested with |
| `experiment1_result.png` | Output plot of Experiment 1 with default settings |

## How to run

```bash
pip install -r requirements.txt
python experiment1_standalone.py
python experiment2_matrix.py
```

Tested on Python 3.12 with NumPy 2.4.4 and Matplotlib 3.10.8.
Each script finishes in under 10 seconds on a laptop.
Experiment 1 writes `experiment1_log.json` and `experiment1_result.png`.
Experiment 2 writes its logs to `exp2_logs/`.

## Example output

Experiment 1 with default settings (k = 1, p = 23, n = 5 to 8).

![Experiment 1 result](experiment1_result.png)

## Verification

The implementation is checked against the invariants listed in Appendix D of the thesis.
Invariant I-1 is reproduced exactly.
For (n, k, m) = (5, 1, 4) the code returns Ω = 288⁴ = 6,879,707,136.

## Known limitation

The default entry points do not reproduce every table in the thesis.
Section 4.1.5 reports v₂(diff) over a grid of n_step and precision values.
Those values are not produced by running the scripts as provided.
The exact configuration behind them has not yet been re-verified.

## Citation

```bibtex
@thesis{yoon2026padic,
  author = {Yoon, Sung-hyuk},
  title  = {On the Origin of Quantum Discreteness in Plasma Electron Scattering:
            A Non-Archimedean Exploration of the Amplituhedron and Its p-adic Structure},
  school = {Dongguk University},
  type   = {Undergraduate thesis},
  year   = {2026}
}
```
