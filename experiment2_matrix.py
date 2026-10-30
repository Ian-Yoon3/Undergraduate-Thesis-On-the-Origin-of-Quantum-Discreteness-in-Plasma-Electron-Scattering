"""
experiment2_matrix.py
=====================
Experiment 2: precision x n_step matrix experiment

purpose:
  (1) effect of precision on depth -> is k=1 depth=64 the theoretical ceiling?
  (2) effect of n_step on depth -> discriminate hypothesis A / B

design:
    precision in {64, 128, 256}
    n_step   in {1, 2, 4}
    k        in {1, 2}
  p=2 fixed, m=4 fixed

discrimination logic:
  - if k=1 depth rises with precision -> k=1's QP=QQ is a precision artifact
  - if k=1 depth stays fixed as precision rises -> structural p-adic completeness of MHV (confirms hypothesis B)
  - if k=2 depth is linear in n_step -> hypothesis A (off-shell deviation)
  - if k=2 depth is independent of n_step -> hypothesis B (structural perturbation)

Run:
  python experiment2_matrix.py
"""

import os
import sys
import time
import json

# reuse the functions of experiment1_standalone.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment1_standalone import (
    run_experiment_1, generate_Z_qq, generate_Z_with_monad,
    is_classical_positive, get_bcfw_cells, compute_total_amplitude,
    compute_padic_error, PadicValuation, FractionMatrix
)


def run_experiment_2(
    k_list=(1, 2),
    precision_list=(64, 128, 256),
    n_step_list=(1, 2, 4),
    n_range_k1=(5, 10),
    n_range_k2=(6, 10),
    p_val=2,
    log_dir="exp2_logs",
):
    """
    Run the precision x n_step matrix.

    Run run_experiment_1 for each (k, precision, n_step) combination.
    Save results in exp2_logs/.
    """
    os.makedirs(log_dir, exist_ok=True)

    all_results = {}
    total_configs = len(k_list) * len(precision_list) * len(n_step_list)
    done = 0
    t_total_start = time.time()

    print("=" * 75)
    print(f"  Experiment 2: precision x n_step matrix")
    print(f"  k={list(k_list)},  precision={list(precision_list)},  n_step={list(n_step_list)}")
    print(f"  {total_configs} combinations total")
    print("=" * 75)

    for k in k_list:
        if k == 1:
            n_min, n_max = n_range_k1
        else:
            n_min, n_max = n_range_k2
        n_list = list(range(n_min, n_max + 1))

        for prec in precision_list:
            for nstep in n_step_list:
                done += 1
                config_key = f"k{k}_prec{prec}_nstep{nstep}"
                log_path = os.path.join(log_dir, f"{config_key}.json")

                print(f"\n[{done}/{total_configs}] {config_key} - n={n_list}")
                print("-" * 70)
                t0 = time.time()

                try:
                    result = run_experiment_1(
                        k=k,
                        n_list=n_list,
                        n_step_monad=nstep,
                        p_val=p_val,
                        precision=prec,
                        log_path=log_path,
                        resume=False,
                    )
                    all_results[config_key] = result
                    elapsed = time.time() - t0
                    print(f"  elapsed: {elapsed:.1f}s")
                except Exception as e:
                    print(f"  failed: {e}")
                    all_results[config_key] = {"error": str(e)}

    # save the full matrix results
    summary = {
        "meta": {
            "k_list": list(k_list),
            "precision_list": list(precision_list),
            "n_step_list": list(n_step_list),
            "p": p_val,
            "total_configs": total_configs,
            "total_elapsed_sec": round(time.time() - t_total_start, 1),
        },
        "results": {}
    }

    for key, result in all_results.items():
        if "error" in result:
            summary["results"][key] = {"error": result["error"]}
            continue
        # per-config summary of (depth, qp_error, omega_qq) by n
        rows_summary = []
        for row in result.get("rows", []):
            rows_summary.append({
                "n":               row.get("n"),
                "depth":           row.get("variation_depth"),
                "qp_error":        row.get("qp_error"),
                "qq_error":        row.get("qq_error"),
                "omega_qq":        row.get("omega_qq"),
                "status":          row.get("status"),
            })
        summary["results"][key] = {
            "k":         result.get("k"),
            "precision": result.get("precision"),
            "n_step":    result.get("n_step"),
            "rows":      rows_summary,
        }

    summary_path = os.path.join(log_dir, "SUMMARY.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, default=str)

    print()
    print("=" * 75)
    print(f"  all done. total elapsed {summary['meta']['total_elapsed_sec']}s")
    print(f"  summary saved: {summary_path}")
    print("=" * 75)

    return summary


def analyze_experiment_2(summary):
    """
    Analyze matrix results -- apply discrimination logic.

    output:
      1. k=1 precision responsiveness (precision artifact vs structural)
      2. k=2 n_step responsiveness (hypothesis A vs B)
      3. conclusion
    """
    print()
    print("=" * 75)
    print("  Analysis 1: precision responsiveness of k=1 depth")
    print("=" * 75)
    print("  n_step=1 fixed; trend of depth(n) as precision varies")
    print()
    print(f"  {'n':>3} | " +
          " | ".join(f"prec={p:>3}" for p in [64, 128, 256]))
    print("  " + "-" * 50)

    # k=1 precision matrix (n_step=1 fixed)
    for n in range(5, 11):
        row_vals = [str(n)]
        for prec in [64, 128, 256]:
            key = f"k1_prec{prec}_nstep1"
            rows = summary["results"].get(key, {}).get("rows", [])
            match = next((r for r in rows if r.get("n") == n), None)
            row_vals.append(f"{match['depth']:>3}" if match and match.get("depth") is not None else "  -")
        print("  " + " | ".join(f"{v:>7}" for v in row_vals))

    # k=1 discrimination
    print()
    print("  discrimination:")
    k1_reactive = False
    for n in [6, 7, 8, 9]:
        depths = []
        for prec in [64, 128, 256]:
            key = f"k1_prec{prec}_nstep1"
            rows = summary["results"].get(key, {}).get("rows", [])
            match = next((r for r in rows if r.get("n") == n), None)
            if match and match.get("depth") is not None:
                depths.append(match["depth"])
        if len(depths) == 3 and (depths[1] > depths[0] or depths[2] > depths[1]):
            k1_reactive = True
            print(f"    n={n}: depth increase confirmed {depths} -> responds to precision")
    if not k1_reactive:
        print("    depth rises in proportion to precision for all n -> precision artifact")
        print("    OR depth stays fixed regardless of precision -> structural completeness (confirms hypothesis B)")

    print()
    print("=" * 75)
    print("  Analysis 2: n_step responsiveness of k=2 depth (discriminate hypothesis A vs B)")
    print("=" * 75)
    print("  precision=128 fixed; trend of depth(n) as n_step varies")
    print()
    print(f"  {'n':>3} | " +
          " | ".join(f"nstep={s:>1}" for s in [1, 2, 4]))
    print("  " + "-" * 45)

    for n in range(6, 11):
        row_vals = [str(n)]
        for nstep in [1, 2, 4]:
            key = f"k2_prec128_nstep{nstep}"
            rows = summary["results"].get(key, {}).get("rows", [])
            match = next((r for r in rows if r.get("n") == n), None)
            row_vals.append(f"{match['depth']:>4}" if match and match.get("depth") is not None else "   -")
        print("  " + " | ".join(f"{v:>7}" for v in row_vals))

    # k=2 n_step response discrimination
    print()
    print("  discrimination:")
    linear_count = 0
    stable_count = 0
    for n in [7, 8, 9, 10]:
        depths = []
        for nstep in [1, 2, 4]:
            key = f"k2_prec128_nstep{nstep}"
            rows = summary["results"].get(key, {}).get("rows", [])
            match = next((r for r in rows if r.get("n") == n), None)
            if match and match.get("depth") is not None:
                depths.append(match["depth"])
        if len(depths) == 3:
            d_stable = len(set(depths)) == 1
            # if depth drops by 1, 2 as n_step doubles 1->2->4, then hypothesis A
            d_linear = (depths[0] - depths[1]) >= 1 and (depths[1] - depths[2]) >= 1
            if d_stable:
                stable_count += 1
                print(f"    n={n}: depth {depths} -- independent of n_step (supports hypothesis B)")
            elif d_linear:
                linear_count += 1
                print(f"    n={n}: depth {depths} -- depth drops as n_step rises (supports hypothesis A)")
            else:
                print(f"    n={n}: depth {depths} -- irregular")

    print()
    print("=" * 75)
    print("  overall conclusion")
    print("=" * 75)

    if stable_count > linear_count:
        print("  * k=2 depth is largely independent of n_step -> supports hypothesis B (structural perturbation)")
    elif linear_count > stable_count:
        print("  * k=2 depth drops as n_step rises -> supports hypothesis A (off-shell deviation)")
    else:
        print("  * k=2 shows both A and B across n ranges -- finer analysis needed")

    print(f"  * hypothesis A evidence: {linear_count} n values")
    print(f"  * hypothesis B evidence: {stable_count} n values")

    return {
        "k1_precision_reactive":    k1_reactive,
        "k2_linear_count":          linear_count,
        "k2_stable_count":          stable_count,
    }


if __name__ == "__main__":
    print("Experiment 2: precision x n_step matrix start\n")

    summary = run_experiment_2(
        k_list         = (1, 2),
        precision_list = (64, 128, 256),
        n_step_list    = (1, 2, 4),
        n_range_k1     = (5, 10),
        n_range_k2     = (6, 10),
        p_val          = 2,
        log_dir        = "exp2_logs",
    )

    # analysis
    verdict = analyze_experiment_2(summary)

    print()
    print("next steps:")
    if verdict["k1_precision_reactive"]:
        print("  * if k=1 also responds to precision, hypothesis B needs review")
    else:
        print("  * k=1 p-adic completeness confirmed -> structural property of the MHV amplitude")

    if verdict["k2_linear_count"] > verdict["k2_stable_count"]:
        print("  * hypothesis A dominant at k=2 -- off-shell deviation is the primary interpretation")
    elif verdict["k2_stable_count"] > verdict["k2_linear_count"]:
        print("  * hypothesis B dominant at k=2 -- structural perturbation is the primary interpretation")
