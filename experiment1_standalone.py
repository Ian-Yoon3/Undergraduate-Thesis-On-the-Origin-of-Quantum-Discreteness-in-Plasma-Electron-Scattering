"""
experiment1_standalone.py
=========================
Amplituhedron Experiment 1: QQ / QP error comparison as particle number n increases

Runs in a standard Python environment (VSCode + WSL Ubuntu) without SageMath.

Required packages:
    pip install numpy matplotlib

Run:
    python experiment1_standalone.py
"""

import gc
import json
import math
import os
import time
from fractions import Fraction
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'  # use English labels


# ================
# 1. Number-system implementation
# ================

class QQField:
    """
    Rational field (exact arithmetic via Fraction).
    Replacement for SageMath's QQ.
    """
    @staticmethod
    def zero():
        return Fraction(0)

    @staticmethod
    def one():
        return Fraction(1)

    @staticmethod
    def from_int(n):
        return Fraction(n)

    @staticmethod
    def to_float(x):
        return float(x)


class PadicValuation:
    """
    p-adic valuation calculator.
    v_p(x) = the number of times p divides x.

    Replacement for SageMath Qp's .valuation().
    For a rational a/b, v_p(a/b) = v_p(a) - v_p(b).
    """
    def __init__(self, p):
        self.p = p

    def valuation(self, x):
        """
        x: Fraction or int
        returns: v_p(x) (integer or inf)
        """
        if x == 0:
            return float('inf')
        if isinstance(x, Fraction):
            return self._val_int(x.numerator) - self._val_int(x.denominator)
        return self._val_int(int(x))

    def _val_int(self, n):
        if n == 0:
            return float('inf')
        n = abs(n)
        count = 0
        while n % self.p == 0:
            n //= self.p
            count += 1
        return count


# =======================
# 2. Matrix operations (exact rational arithmetic)
# =======================

class FractionMatrix:
    """
    Exact matrix operations based on Fraction.
    Replacement for SageMath matrix(QQ, ...).
    Computes the determinant as an exact rational (no floating-point error).
    """
    def __init__(self, rows):
        """
        rows: list of list of Fraction (int/float are converted to Fraction)
        """
        self.data = [[Fraction(x) for x in row] for row in rows]
        self.nrows = len(rows)
        self.ncols = len(rows[0]) if rows else 0

    def matrix_from_rows_and_columns(self, row_indices, col_indices):
        """Extract a submatrix from the given rows/columns."""
        rows = []
        for r in row_indices:
            rows.append([self.data[r][c] for c in col_indices])
        return FractionMatrix(rows)

    def matrix_from_columns(self, col_indices):
        """Extract a submatrix from the given columns."""
        return self.matrix_from_rows_and_columns(range(self.nrows), col_indices)

    def det(self):
        """
        Fraction determinant via LU decomposition.
        Guarantees an exact rational result.
        """
        n = self.nrows
        if n != self.ncols:
            raise ValueError(
                f"det() applies only to square matrices. size: {self.nrows}x{self.ncols}"
            )
        if n == 0:
            return Fraction(1)
        if n == 1:
            return self.data[0][0]
        if n == 2:
            return (self.data[0][0] * self.data[1][1]
                    - self.data[0][1] * self.data[1][0])

        # Gaussian elimination with Fraction (exact)
        M = [row[:] for row in self.data]
        det = Fraction(1)
        for col in range(n):
            # pivot selection
            pivot_row = None
            for row in range(col, n):
                if M[row][col] != 0:
                    pivot_row = row
                    break
            if pivot_row is None:
                return Fraction(0)
            if pivot_row != col:
                M[col], M[pivot_row] = M[pivot_row], M[col]
                det = -det

            det *= M[col][col]
            pivot = M[col][col]
            for row in range(col + 1, n):
                if M[row][col] == 0:
                    continue
                factor = M[row][col] / pivot
                for c in range(col, n):
                    M[row][c] -= factor * M[col][c]

        return det


# ========================
# 3. Momentum-twistor Z generation (moment curve)
# =======================
def generate_Z_qq(n, k, m=4):
    """
    Rational track: Z_a = (1, a, a^2, ..., a^(k+m-1))
    returns: FractionMatrix (k+m) x n
    """
    rows = []
    for a in range(1, n + 1):
        row = [Fraction(a ** j) for j in range(k + m)]
        rows.append(row)
    # transpose to (k+m) x n
    transposed = [[rows[a][j] for a in range(n)] for j in range(k + m)]
    return FractionMatrix(transposed)


# =======================
# 4. Positivity condition check
# ======================
def is_classical_positive(Z_matrix, k, m=4):
    """
    Check ordered k x k minor > 0 -- definition eq.(2.4) of the paper.

    [Theoretical basis]
    The positivity condition of the amplituhedron requires that only
    'ordered' minors with a1 < a2 < ... < ak be positive.
    Postnikov (2006)'s sufficient condition is likewise based on ordered minors.

    [Cause of a previous bug]
    Using cyclic-consecutive indices (including wrap-around)
    includes reversed combinations such as cols=(n-1, 0), causing det < 0.
    The moment curve guarantees positivity only for ordered minors;
    the sign of a wrap-around minor is arbitrary depending on k, n.

    [Fix]
    iterate a1 over range(n - k + 1) -> always guarantees a1 < a2 < ... < ak.
    no wrap-around -> passes correctly for all k=1,2,3,...

    [Squareness]
    use the top k rows x k columns for the k x k minor.
    """
    n_cols = Z_matrix.ncols

    if k == 0:
        return True  # G(0,n): always positive

    # ordered: a1 < a2 < ... < ak, no wrap
    for a1 in range(n_cols - k + 1):
        cols = list(range(a1, a1 + k))
        # top k rows x k cols -> guarantees k x k square
        sub = Z_matrix.matrix_from_rows_and_columns(list(range(k)), cols)
        val = sub.det()
        if val <= 0:
            return False
    return True


# ======================
# 5. BCFW cell decomposition (with memoization)
# ======================

def get_bcfw_cells(n, k, m=4, current_labels=None, cache=None):
    """
    BCFW recursive cell decomposition.
    Memoization: cached with (labels_tuple, k) as key.
    """
    if cache is None:
        cache = {}
    if current_labels is None:
        current_labels = list(range(n))

    cache_key = (tuple(current_labels), k)
    if cache_key in cache:
        return cache[cache_key]

    # physical validity
    if n < k + m:
        raise ValueError(f"physically invalid: n={n} < k+m={k+m}")

    # base case 1: k == 0 -> MHV
    if k == 0:
        result = {"type": "MHV", "indices": tuple(current_labels), "k": 0}
        cache[cache_key] = result
        return result

    # base case 2: n == k + m -> simplex
    if n == k + m:
        result = {
            "type":    "simplex",
            "indices": tuple(sorted(current_labels)),
            "k":       k
        }
        cache[cache_key] = result
        return result

    all_cells = []
    for i in range(2, n - 1):
        k_L_min = max(0, k - (n - i + 1))
        k_L_max = min(k - 1, i - 2)

        for k_L in range(k_L_min, k_L_max + 1):
            k_R = k - 1 - k_L
            l_labels = current_labels[:i + 1]
            r_labels = [current_labels[0]] + current_labels[i:]

            if len(l_labels) < k_L + m:
                continue
            if len(r_labels) < k_R + m:
                continue

            l_res = get_bcfw_cells(len(l_labels), k_L, m, l_labels, cache)
            r_res = get_bcfw_cells(len(r_labels), k_R, m, r_labels, cache)

            l_list = l_res if isinstance(l_res, list) else [l_res]
            r_list = r_res if isinstance(r_res, list) else [r_res]

            for l_cell in l_list:
                for r_cell in r_list:
                    all_cells.append({
                        "type":    "bcfw",
                        "L_part":  l_cell,
                        "R_part":  r_cell,
                        "split_at": i,
                        "k_map":   (k_L, k_R),
                        "k":       k
                    })

    cache[cache_key] = all_cells
    return all_cells


# =======================
# 6. Canonical form computation
# =======================
def _slice_Z(Z_matrix, k_sub, indices, m=4):
    """
    slice Z rows to k_sub -> guarantees (k_sub+m) x (k_sub+m) square.
    """
    k_plus_m = k_sub + m
    n_rows = Z_matrix.nrows

    if n_rows == k_plus_m:
        return Z_matrix.matrix_from_columns(indices)
    elif n_rows > k_plus_m:
        return Z_matrix.matrix_from_rows_and_columns(
            list(range(k_plus_m)), indices
        )
    else:
        raise ValueError(
            f"Z row count({n_rows}) < k_sub+m({k_plus_m}). "
            f"k_sub={k_sub}, indices={indices}"
        )


def compute_simplex_canonical_form(cell, Z_matrix, m=4):
    """
    Compute the canonical form of a simplex cell.
    Omega = det(Z_I)^m / prod(cyclic consecutive k_sub x k_sub minors)

    [Squareness note]
    - numerator det: (k_sub+m)x(k_sub+m) -> use _slice_Z
    - denominator minor: k_sub x k_sub -> use only the top k_sub rows
    """
    indices = list(cell["indices"])
    k_sub   = cell["k"]
    n_local = len(indices)  # == k_sub + m

    # numerator: det(Z_sq)^m -- (k_sub+m)x(k_sub+m) square
    Z_sq    = _slice_Z(Z_matrix, k_sub, indices, m)
    det_val = Z_sq.det()
    numerator = det_val ** m

    # denominator: product of cyclic-consecutive k_sub x k_sub minors
    # must be k_sub rows x k_sub cols -> k_sub x k_sub square
    denominator = Fraction(1)
    for start in range(n_local):
        local_cols  = [(start + offset) % n_local for offset in range(k_sub)]
        global_cols = [indices[c] for c in local_cols]

        # top k_sub rows x k_sub cols -> guarantees k_sub x k_sub square
        sub       = Z_matrix.matrix_from_rows_and_columns(
                        list(range(k_sub)), global_cols
                    )
        minor_val = sub.det()
        if minor_val == 0:
            return Fraction(0)
        denominator *= minor_val

    if denominator == 0:
        return Fraction(0)
    return numerator / denominator


def compute_cell_canonical_form(cell, Z_matrix, m=4):
    """
    Recursively compute the canonical form by cell type.
    returns: Fraction
    """
    cell_type = cell.get("type")

    if cell_type == "MHV":
        return Fraction(1)

    if cell_type == "simplex":
        return compute_simplex_canonical_form(cell, Z_matrix, m)

    if cell_type == "bcfw":
        total = Fraction(0)
        l_cells = cell["L_part"] if isinstance(cell["L_part"], list) else [cell["L_part"]]
        r_cells = cell["R_part"] if isinstance(cell["R_part"], list) else [cell["R_part"]]
        for l_cell in l_cells:
            for r_cell in r_cells:
                l_val = compute_cell_canonical_form(l_cell, Z_matrix, m)
                r_val = compute_cell_canonical_form(r_cell, Z_matrix, m)
                total += l_val * r_val
        return total

    return Fraction(0)


def compute_total_amplitude(bcfw_cells, Z_matrix, m=4):
    """
    Sum the canonical forms of all BCFW cells -> total amplitude Omega.
    returns: Fraction
    """
    cells = bcfw_cells if isinstance(bcfw_cells, list) else [bcfw_cells]
    total = Fraction(0)
    for cell in cells:
        total += compute_cell_canonical_form(cell, Z_matrix, m)
    return total


# =======================
# 7. p-adic error computation
# ======================

def compute_padic_error(omega_qq, omega_qp, p, prec):
    """
    Compute the p-adic error via the relative-depth method.

    [Key idea]
    with eps = Fraction(n_step, p^N), v_p(eps) = -N.
    diff = omega_qp - omega_qq ~ (d omega/d eps) * eps

    raw depth = v_p(diff) = v_p(d omega/d eps) + v_p(eps)
                            = v_p(partial omega/partial eps) - N
    -> depends on N, so not comparable across n.

    relative depth = v_p(diff) - v_p(eps)
               = v_p(diff) + N   (since v_p(eps) = -N)
               = v_p(partial omega/partial eps)
    -> p-adic structure of the amplitude itself, independent of N (the eps setting)
    -> a meaningful quantity that varies with n

    QP_error = p^(-relative_depth)
    the larger relative_depth is, the smaller QP_error is
    (stronger p-adic structure of the amplitude -> perturbation is p-adically aligned)

    Parameters
    ----------
    omega_qq : Fraction  rational amplitude
    omega_qp : Fraction  perturbed amplitude (includes eps = 1/p^N)
    p        : int       prime
    prec     : int       precision (N = prec//2)

    Returns
    -------
    (relative_depth, qp_error, raw_depth)
    """
    pval = PadicValuation(p)
    N    = prec // 2           # v_p(eps) = -N

    diff = omega_qp - omega_qq

    if diff == 0:
        # exact match: use the theoretical ceiling prec as the relative depth
        return prec, 0.0, prec

    raw_depth = pval.valuation(diff)

    if raw_depth == float('inf'):
        return prec, 0.0, prec

    raw_int = int(raw_depth)

    # relative depth = v_p(diff) - v_p(eps) = raw_depth - (-N) = raw_depth + N
    relative_depth = raw_int + N

    if relative_depth <= 0:
        # perturbation has no p-adic structure -> QP_error = 1 (maximum error)
        qp_error = 1.0
    else:
        qp_error = float(p) ** (-relative_depth)

    return relative_depth, qp_error, raw_int


def generate_Z_with_monad(n, k, m=4, n_step=1, p=23, prec=1024):
    """
    Z matrix with monad perturbation applied (QP track).

    [Bug history]
    v1: eps = p^(prec//2) (large integer) -> Z information is buried in eps
    v2: eps = 1/p^N (small rational, uniform) -> fraction reduction fixes v_p(diff) = 0

    [Root cause of v2]
    adding the same eps to all Z entries means that, when computing omega_qp - omega_qq,
    fraction reduction cancels 23^N from both numerator and denominator.
    result: v_p(diff) = 0 -> relative_depth = 128 (a constant independent of n)

    [v3 fix: per-entry distinct eps]
    eps_{row, col} = eps_base \times (row+1) \times (col+1)
    = (n_step / p^N) \times (row+1) \times (col+1)

    reason:
    - (row+1), (col+1) are integers in 1..(k+m), 1..n
    - since n < p=23, (col+1) <= n < 23 -> no multiple of 23
    - (row+1) <= k+m <= 7 -> no multiple of 23
    - hence v_p(eps_{row,col}) = v_p(eps_base) = -N
    - but since each entry's eps differs, in the determinant expansion
      the eps factors of each term differ, so cancellation is incomplete
    - -> v_p of diff = omega_qp - omega_qq acquires n-dependence

    physical meaning:
    (row+1): weight along the moment-curve order (higher orders perturbed more)
    (col+1): weight along the particle number (later particles perturbed more)
    this models a non-uniform off-shell deviation as particles are added.
    """
    eps_base = Fraction(n_step, p ** (prec // 2))

    rows = []
    for a in range(1, n + 1):          # a = col+1 (1-indexed)
        row = []
        for j in range(k + m):         # j = row index (0-indexed)
            eps_ij = eps_base * (j + 1) * a   # per-entry distinct eps
            row.append(Fraction(a ** j) + eps_ij)
        rows.append(row)
    transposed = [[rows[a][j] for a in range(n)] for j in range(k + m)]
    return FractionMatrix(transposed)


# ======================
# 8. Experiment 1 main function
# ======================
def run_experiment_1(
    k=1,
    n_list=None,
    n_step_monad=1,
    p_val=23,
    precision=1024,
    log_path="experiment1_log.json",
    resume=True,
    m=4
):
    """
    Experiment 1: compare QQ / QP error as particle number n increases.

    Parameters
    ----------
    k            : NMHV degree (default 1)
    n_list       : list of n to test. If None, defaults to [k+m, ..., k+m+5].
    n_step_monad : monad perturbation strength
    p_val        : p-adic prime (default 23, fixed throughout)
    precision    : upper bound for variation_depth (default 1024)
    log_path     : path to the JSON log file
    resume       : if True, resume from an existing log
    m            : spacetime dimension (fixed 4)
    """
    n_min_valid = k + m
    if n_list is None:
        n_list = list(range(n_min_valid, n_min_valid + 6))

    invalid = [n for n in n_list if n < n_min_valid]
    if invalid:
        raise ValueError(f"invalid n included: {invalid}. minimum={n_min_valid}")

    # -- load existing log --
    completed = {}
    omega_ref = None

    if resume and os.path.exists(log_path):
        try:
            with open(log_path, "r", encoding="utf-8") as f:
                log_data = json.load(f)
            for row in log_data.get("rows", []):
                completed[row["n"]] = row
            omega_ref = log_data.get("omega_ref")
            print(f"[Resume] '{log_path}' -> {len(completed)} results loaded.")
            remaining = [n for n in n_list if n not in completed]
            print(f"  Done: {sorted(completed.keys())}  Remaining: {remaining}")
            if not remaining:
                print("  All done. Returning cached results.")
                rows = [completed[n] for n in sorted(completed) if n in n_list]
                return _build_result(k, p_val, precision, n_step_monad,
                                     n_list, rows, omega_ref, log_path)
            print()
        except Exception as e:
            print(f"[Warning] Log read failed ({e}). Starting fresh.")
            completed = {}
            omega_ref = None

    # -- print experiment header --
    print("=" * 70)
    print("  Experiment 1: QQ vs QP Error Rate by Particle Number n")
    print(f"  k={k}, m={m}, p={p_val} (fixed), precision={precision} (fixed)")
    print(f"  n_list={n_list}, n_step_monad={n_step_monad}")
    print(f"  Log: {log_path}")
    print("=" * 70)
    print(f"  {'n':>3} | {'cells':>5} | {'Omega_QQ':>18} | "
          f"{'QQ_err':>12} | {'depth':>6} | {'QP_err':>14} | {'time':>6} | status")
    print("-" * 100)

    rows_all = list(completed.values())
    n_min_in_list = min(n_list)

    for n in n_list:
        if n in completed:
            r = completed[n]
            print(f"  {n:>3} | [skipped - already in log]")
            continue

        row = {
            "n":         n,
            "p_used":    p_val,
            "prec_used": precision,
            "status":    "ok",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        t0 = time.time()
        bcfw_cells = None
        Z_qq       = None
        Z_qp       = None

        try:
            # -- 1. generate positive Z matrix --
            Z_qq = generate_Z_qq(n, k, m)

            # verify positivity condition
            if not is_classical_positive(Z_qq, k):
                raise ValueError("Z_qq does not satisfy the positivity condition.")

            # -- 2. monad-perturbed Z (QP track) -- pass precision to fix eps direction
            Z_qp = generate_Z_with_monad(n, k, m, n_step_monad, p_val, precision)

            # -- 3. BCFW cell decomposition --
            cache = {}
            bcfw_cells = get_bcfw_cells(n, k, m, cache=cache)
            n_cells = len(bcfw_cells) if isinstance(bcfw_cells, list) else 1
            row["n_cells"] = n_cells
            del cache  # release cache immediately

            # -- 4. amplitude computation (QQ track) --
            omega_qq_frac = compute_total_amplitude(bcfw_cells, Z_qq, m)
            omega_qq_f    = float(omega_qq_frac)
            row["omega_qq"] = omega_qq_f

            # -- 5. amplitude computation (QP track) --
            omega_qp_frac = compute_total_amplitude(bcfw_cells, Z_qp, m)
            omega_qp_f    = float(omega_qp_frac)
            row["omega_qp_real"] = omega_qp_f

            # -- 6. p-adic variation depth (relative-depth method) --
            rel_depth, qp_err, raw_depth = compute_padic_error(
                omega_qq_frac, omega_qp_frac, p_val, precision
            )
            row["variation_depth"] = rel_depth   # relative depth (captures n-dependence)
            row["raw_depth"]       = raw_depth   # raw depth (for reference)
            row["qp_error"]        = qp_err

            # -- 7. QQ error computation --
            # reference: Omega_QQ of the first successful n (proceed to next n if n_min fails)
            if omega_ref is None:
                omega_ref       = omega_qq_f if omega_qq_f != 0 else 1.0
                row["qq_error"] = 0.0
            elif omega_ref != 0:
                row["qq_error"] = abs(omega_qq_f - omega_ref) / abs(omega_ref)
            else:
                row["qq_error"] = float('inf')

        except Exception as e:
            row["status"] = f"ERROR: {e}"
            for key in ["n_cells", "omega_qq", "omega_qp_real",
                        "variation_depth", "qq_error", "qp_error"]:
                row.setdefault(key, None)

        # -- 8. release memory immediately --
        del bcfw_cells
        del Z_qq
        del Z_qp
        gc.collect()

        row["elapsed_sec"] = round(time.time() - t0, 2)

        # -- 9. one-line output --
        omega_s  = f"{row.get('omega_qq'):.6e}"    if row.get('omega_qq')       is not None else "N/A"
        qq_s     = f"{row.get('qq_error'):.4e}"    if row.get('qq_error')       is not None else "N/A"
        qp_s     = f"{row.get('qp_error'):.4e}"    if row.get('qp_error')       is not None else "N/A"
        dp_s     = str(row.get('variation_depth',  'N/A'))
        cells_s  = str(row.get('n_cells', '-'))

        print(f"  {n:>3} | {cells_s:>5} | {omega_s:>18} | "
              f"{qq_s:>12} | {dp_s:>6} | {qp_s:>14} | "
              f"{row['elapsed_sec']:>5}s | {row['status']}")

        # -- 10. save log immediately --
        rows_all.append(row)
        completed[n] = row
        _save_log(log_path, k, p_val, precision, n_step_monad,
                  n_list, rows_all, omega_ref)

    rows_sorted = [completed[n] for n in sorted(completed) if n in n_list]
    return _build_result(k, p_val, precision, n_step_monad,
                         n_list, rows_sorted, omega_ref, log_path)


# ==============================
# 9. log-saving / result-assembly helpers
# ==============================
def _save_log(log_path, k, p, precision, n_step, n_list, rows, omega_ref):
    """Save JSON immediately after each n completes."""
    safe_rows = []
    for row in rows:
        safe_row = {}
        for key, val in row.items():
            if val == float('inf'):
                safe_row[key] = "Inf"
            elif val != val:  # NaN
                safe_row[key] = "NaN"
            else:
                try:
                    json.dumps(val)
                    safe_row[key] = val
                except Exception:
                    safe_row[key] = str(val)
        safe_rows.append(safe_row)

    data = {
        "experiment": "experiment_1",
        "k":          k,
        "p":          p,
        "precision":  precision,
        "n_step":     n_step,
        "n_list":     n_list,
        "omega_ref":  omega_ref,
        "saved_at":   time.strftime("%Y-%m-%d %H:%M:%S"),
        "rows":       safe_rows,
    }
    try:
        with open(log_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"  [Log save failed] {e}")


def _build_result(k, p, precision, n_step, n_list, rows, omega_ref, log_path):
    return {
        "k":         k,
        "p":         p,
        "precision": precision,
        "n_step":    n_step,
        "n_list":    n_list,
        "omega_ref": omega_ref,
        "log_path":  log_path,
        "rows":      rows,
    }


# ===========================
# 10. result output
# ===========================

def print_experiment_1(results):
    """Print a summary of Experiment 1 results."""
    k    = results["k"]
    p    = results["p"]
    prec = results["precision"]

    print()
    print("=" * 80)
    print("  Experiment 1 Results: QQ vs QP Error Rate")
    print(f"  k={k}, p={p} (fixed), precision={prec} (fixed), m=4")
    print(f"  Reference Omega_QQ (n_min) = "
          f"{results['omega_ref']:.6g}" if results.get('omega_ref') is not None
          else "  Reference Omega_QQ (n_min) = N/A")
    print("=" * 80)
    print(f"  {'n':>3} | {'cells':>5} | {'Omega_QQ':>16} | "
          f"{'QQ_error':>12} | {'depth':>6} | {'QP_error':>14} | interpretation")
    print("-" * 90)

    for row in results["rows"]:
        n = row["n"]
        if row.get("status", "ok") != "ok":
            print(f"  {n:>3} | ERROR: {row['status']}")
            continue

        omega_s = f"{row.get('omega_qq', 0):.6g}"
        qq_s    = f"{row.get('qq_error', 0):.4e}"   if row.get('qq_error')  is not None else "N/A"
        qp_s    = f"{row.get('qp_error', 0):.4e}"   if row.get('qp_error')  is not None else "N/A"
        dp_s    = str(row.get('variation_depth', 'N/A'))

        # interpretation
        qe = row.get('qq_error')
        pe = row.get('qp_error')
        if pe is not None and qe is not None:
            if pe == 0.0:
                interp = "QP == QQ (exact)"
            elif pe < 1e-10:
                interp = "QP near-exact"
            elif pe < qe:
                interp = "QP < QQ error"
            elif pe > qe:
                interp = "QP > QQ error"
            else:
                interp = "QP = QQ error"
        else:
            interp = ""

        print(f"  {n:>3} | {row.get('n_cells','-'):>5} | {omega_s:>16} | "
              f"{qq_s:>12} | {dp_s:>6} | {qp_s:>14} | {interp}")

    print("=" * 80)

    # trend summary
    valid = [r for r in results["rows"]
             if r.get("status", "ok") == "ok"
             and r.get("qq_error") is not None
             and r.get("qp_error") is not None]
    if len(valid) < 2:
        return

    qq_errs = [r["qq_error"] for r in valid]
    qp_errs = [r["qp_error"] for r in valid]
    depths  = [r["variation_depth"] for r in valid if r.get("variation_depth") is not None]

    def trend(vals):
        if all(vals[i] <= vals[i+1] for i in range(len(vals)-1)):
            return "monotone increasing"
        elif all(vals[i] >= vals[i+1] for i in range(len(vals)-1)):
            return "monotone decreasing"
        return "irregular"

    print()
    print("  [Summary]")
    print(f"  QQ_error trend : {trend(qq_errs)}")
    print(f"  QP_error trend : {trend(qp_errs)}")
    qp_lt = sum(1 for r in valid if r["qp_error"] < r["qq_error"])
    print(f"  QP < QQ error  : {qp_lt}/{len(valid)} cases")
    if depths:
        print(f"  depth range    : {min(depths)} ~ {max(depths)}")
        if len(set(depths)) == 1:
            print(f"  depth stable ({depths[0]}) -> structural perturbation (Hypothesis B)")


# ===========================
# 11. plot output
# ============================

def plot_experiment_1(results, save_path=None):
    """
    Grouped bar chart + log y-axis.
    - blue bars: QQ_error
    - orange bars: QP_error
    - green line (right axis): variation_depth
    """
    valid = [r for r in results["rows"]
             if r.get("status", "ok") == "ok"
             and r.get("qq_error") is not None
             and r.get("qp_error") is not None]

    if not valid:
        print("No valid data for plot.")
        return

    ns        = [r["n"]               for r in valid]
    qq_errors = [max(r["qq_error"],  1e-300) for r in valid]
    qp_errors = [max(r["qp_error"],  1e-300) for r in valid]
    depths    = [r.get("variation_depth") for r in valid]

    fig, ax1 = plt.subplots(figsize=(max(8, len(ns) * 1.4), 5.5))

    x = np.arange(len(ns))
    w = 0.30

    bars_qq = ax1.bar(x - w/2, qq_errors, w,
                      label="QQ_error (rational)", color="#3B82F6", alpha=0.85)
    bars_qp = ax1.bar(x + w/2, qp_errors, w,
                      label=f"QP_error (p={results['p']}-adic)", color="#F97316", alpha=0.85)

    ax1.set_yscale("log")
    ax1.set_xlabel("Particle number n", fontsize=11)
    ax1.set_ylabel("Error rate (log scale)", fontsize=11)
    ax1.set_xticks(x)
    ax1.set_xticklabels([str(n) for n in ns])
    ax1.grid(axis="y", alpha=0.3, linestyle="--")
    ax1.set_ylim(bottom=1e-300 * 0.1)

    # value labels
    for bar, val in zip(bars_qq, [r["qq_error"] for r in valid]):
        if val is not None and val > 1e-100:
            ax1.text(bar.get_x() + bar.get_width()/2,
                     bar.get_height() * 2,
                     f"{val:.1e}", ha="center", va="bottom",
                     fontsize=7, color="#1E40AF")
    for bar, val in zip(bars_qp, [r["qp_error"] for r in valid]):
        if val is not None and val > 1e-100:
            ax1.text(bar.get_x() + bar.get_width()/2,
                     bar.get_height() * 2,
                     f"{val:.1e}", ha="center", va="bottom",
                     fontsize=7, color="#C2410C")

    # secondary axis: variation_depth
    ax2 = ax1.twinx()
    valid_d = [(i, d) for i, d in enumerate(depths) if d is not None]
    if valid_d:
        xi, di = zip(*valid_d)
        ax2.plot(list(xi), list(di), "g-o", linewidth=1.8,
                 markersize=5, label="variation_depth", alpha=0.85)
    ax2.set_ylabel("variation_depth  v_p(Omega_QP - Omega_QQ)",
                   fontsize=10, color="green")
    ax2.tick_params(axis="y", labelcolor="green")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=9)

    title = (f"Experiment 1: QQ vs QP Error Rate by n\n"
             f"k={results['k']}, m=4, p={results['p']} (fixed), "
             f"precision={results['precision']}, "
             f"n_step={results['n_step']}")
    ax1.set_title(title, fontsize=10, pad=10)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Graph saved: {save_path}")
    else:
        plt.savefig("experiment1_result.png", dpi=150, bbox_inches="tight")
        print("Graph saved: experiment1_result.png")
    plt.show()
    plt.close()


# =========================
# 12. estimate n limit
# =========================

def estimate_n_limit(k=1, n_max=15, m=4):
    """Estimate the feasible n limit by counting BCFW cells only."""
    cache = {}

    def count(n, k):
        key = (n, k)
        if key in cache: return cache[key]
        if n < k + m: return 0
        if k == 0: return 1
        if n == k + m: return 1
        total = 0
        for i in range(2, n - 1):
            k_L_min = max(0, k - (n - i + 1))
            k_L_max = min(k - 1, i - 2)
            for k_L in range(k_L_min, k_L_max + 1):
                k_R = k - 1 - k_L
                if i+1 >= k_L+m and n-i+1 >= k_R+m:
                    total += count(i+1, k_L) * count(n-i+1, k_R)
        cache[key] = total
        return total

    print(f"\n  n limit estimate (k={k}, m={m})")
    print(f"  {'n':>4} | {'cells':>8} | {'verdict':>20}")
    print("  " + "-" * 40)
    prev = 1
    for n in range(k+m, n_max+1):
        c = count(n, k)
        g = c / prev if prev > 0 else 0
        prev = max(c, 1)
        v = ("recommended" if c <= 50 else
             "feasible"    if c <= 500 else
             "slow"        if c <= 5000 else "too slow")
        print(f"  {n:>4} | {c:>8,} | {v:>20}  (x{g:.1f})")


# =====================
# 13. entry point
# =======================

if __name__ == "__main__":
    print()
    print("Amplituhedron Experiment 1 - Standalone Python")
    print("p=23 fixed, comparing QQ (rational) vs QP (p-adic) error rates")
    print()

    # -- preliminary check: n limit --
    estimate_n_limit(k=1, n_max=12)
    print()

    # -- run experiment --
    # specify n_list directly; run only the desired n values.
    # precision: 1024 (default) or 2048 (more precise, slower)
    results = run_experiment_1(
        k            = 1,
        n_list       = [5, 6, 7, 8],   # <- specify the desired n list directly
        n_step_monad = 1,
        p_val        = 23,
        precision    = 1024,
        log_path     = "experiment1_log.json",
        resume       = True,
    )

    # -- print results --
    print_experiment_1(results)

    # -- save plot --
    plot_experiment_1(results, save_path="experiment1_result.png")


