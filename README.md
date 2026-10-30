# p-adic Structure of the Amplituhedron

Source code for the undergraduate thesis
**On the Origin of Quantum Discreteness in Plasma Electron Scattering. A Non-Archimedean Exploration of the Amplituhedron and Its p-adic Structure**
(Sung-hyuk Yoon, Dongguk University, 2026. Supervised by Dr. Jeon ByungHoon.)

📄 [Read the full thesis (PDF, 106 pages)](thesis.pdf)

The thesis is also featured on the author's [LinkedIn profile](https://www.linkedin.com/in/sunghyuk-yoon-3a1319383/).

## Motivation

As semiconductor processes approach the nanometer scale, stochastic fluctuation sets a limit that precise control does not remove.
This thesis asks a question one step behind that limit.
Probability in the microscopic world may not be a feature of nature itself.
It may instead be noise from describing nature in the language of real numbers.
The amplituhedron, a geometric object that encodes scattering amplitudes, gives a setting where this idea can be tested by exact computation.

## Overview

The code computes the canonical form of the amplituhedron in two parallel tracks.
The QQ track works in exact rational arithmetic.
The QP track adds a p-adically scaled perturbation to the same input.
Both tracks are exact, so the comparison is not a comparison of floating-point error.
The difference between the two tracks is read through its p-adic valuation instead of its magnitude.

![Simulation flowchart](simulation_algorithm.png)

*Flowchart of the simulation (left) and its correspondence with the amplituhedron construction of Arkani-Hamed and Trnka (right).*

## Key findings

The QP track does not compute the amplitude more precisely than the QQ track.
Both tracks carry rational entries, and neither approximates anything.

Read through the p-adic valuation, the two tracks show a definite structure.
At k = 1 the amplitude reduces to the number of cells for n ≥ 6, so the two tracks return the same value.
At k = 2 the behavior divides sharply between n = 7 and n = 8.
Below this threshold the valuation of the difference equals the valuation of the amplitude.

| (n, k) | v₂(Ω_QP − Ω_QQ) | v₂(Ω_QQ) |
|---|---|---|
| (5, 1) | 20 | 20 |
| (6, 2) | 32 | 32 |
| (7, 2) | 21 | 21 |

The agreement does not depend on the size of the perturbation.
It holds at a single value across thirty combinations of n_step and precision.
It disappears when the perturbation is scaled to a prime other than the one used for the valuation.

The structure has a traceable origin.
On the moment curve, the numerator of each cell is a Vandermonde determinant.
The prime content of the amplitude is inherited from the differences of small integers in that determinant.

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

![BCFW recursion tree for (7, 2)](bcfwtree.png)

*BCFW recursion tree for (n, k) = (7, 2). Each leaf is a base case, either an MHV cell or a simplex cell.*

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
| `simulation_algorithm.png`, `bcfwtree.png`, `experiment1_result.png` | Flowchart, recursion tree and output plot |
| `thesis.pdf`, `summary.pdf` | Full text of the thesis and a five-page summary |

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

## Reading the output

The scripts report a relative depth, not the raw valuation.
The relative depth equals v_p(Ω_QP − Ω_QQ) + prec/2.
Subtracting prec/2 recovers the values in the table above.
For example, Experiment 2 reports a depth of 52, 84 and 148 for (5, 1) at prec = 64, 128 and 256.
Each gives v₂ = 20 after the subtraction.
When the two tracks agree exactly, the scripts report the depth as prec by convention.

Experiment 1 with default settings (k = 1, p = 23, n = 5 to 8) produces the plot below.

![Experiment 1 result](experiment1_result.png)

## Verification

The implementation is checked against the invariants listed in Appendix D of the thesis.
Invariant I-1 is reproduced exactly.
For (n, k, m) = (5, 1, 4) the code returns Ω = 288⁴ = 6,879,707,136.

## Thesis summary

A five-page summary of the thesis is available for a quicker read.

📄 [Read the summary (PDF)](summary.pdf)

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
