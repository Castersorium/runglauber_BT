"""
pA_from_tglauber_uproot.py

Requirements:
  pip install uproot numpy matplotlib scipy

Usage:
  python pA_from_tglauber_uproot.py /path/to/your.root

说明:
- 自动检查树 nt_p_Au 的分支名并尝试以下检测顺序来取得 "每个投射核子 (projectile-side) 的碰撞次数 n_i":
    1) 优先查找明显的 per-nucleon / jagged-array 分支，如 'NcollA', 'NcollB', 'fNcollA_arr', 'collisionsA', 'ncollA' 等。
    2) 如果没有 per-nucleon arrays：对 p+A 事件（NpartA≈1）使用 event-level scalar Ncoll/NcollA/NcollA_scalar 作为该 projectile 的 ν。
    3) 再没有的话，将尝试用 Ncoll（总二体碰撞数）作为近似（仅适用于 p+A 且当 event NpartA==1 时）。
- 采样采用 exponential 分布: Delta y_i ~ Exp(scale = 1/alpha)
  允许第 1 次与 subsequent 有不同 alpha (alpha1, alpha_rest).
- 输出:
  - 一个 dN/dy 直方图 (projectile-side final rapidity)
  - 若需要，可扩展为拟合 experimental data 的接口
"""

import sys
import os
import numpy as np
import uproot
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.optimize import minimize

# --------------------------- 用户配置 ---------------------------
# 输入 ROOT 文件（从命令行传入，第一个参数）
if len(sys.argv) < 2:
    print("Usage: python pA_from_tglauber_uproot.py /path/to/file.root [options]")
    sys.exit(1)

FILENAME = sys.argv[1]
TREENAME = "nt_p_Au"   # 你给出的 tree 名称
# beam rapidity (调整为你实验能量对应的数值)
YBEAM = 5.36

# alpha 参数初始值（可用于拟合）
ALPHA1 = 1.0       # first collision alpha -> mean loss 1/alpha1 ~ 1.0
ALPHAREST = 3.0    # subsequent collisions alpha -> mean loss 1/alpha_rest ~ 1/3

# histogram settings
NBINS = 200
YRANGE = (-5.0, 5.0)

# centrality selection example (optional) - 使用 Npart 或 B 来筛选
CENTRALITY_CUT = None   # e.g. ( "Npart", 2, 10 )  -> 保留 Npart in [2,10]
# ------------------------ end user config ------------------------

def list_tree_branches(filename, treename):
    with uproot.open(filename) as f:
        if treename not in f:
            raise KeyError(f"Tree {treename} not found in {filename}. Available keys: {list(f.keys())}")
        tree = f[treename]
        print("Branches in tree (name : interpretation if available):")
        for k, v in tree.items():
            print(f"  {k}  -> {v}")
        return tree

def detect_candidate_branches(tree):
    """
    返回可能包含 per-nucleon 或 event-level Ncoll 的分支名集合
    """
    names = list(tree.keys())
    names_lower = [n.lower() for n in names]

    # 常见候选模式（按优先级）
    per_nucleon_patterns = [
        "ncoll_a", "ncollb", "ncoll_a_arr", "ncollb_arr", "ncoll_a[]", "fNcollA", "fNcollB",
        "collisions_a", "collisionsb", "collisions_a_arr", "collisions_b_arr",
        "collisions_a[]", "collisions", "ncoll", "ncolls"
    ]
    event_level_candidates = [
        "ncoll", "npart", "nparta", "npartb", "ncollpp", "ncollpn", "ncollnn"
    ]

    found_per = []
    found_event = []
    for i,n in enumerate(names_lower):
        for p in per_nucleon_patterns:
            if p in n:
                found_per.append(names[i])
        for p in event_level_candidates:
            if p in n and names[i] not in found_per:
                found_event.append(names[i])

    # make unique in order
    found_per = list(dict.fromkeys(found_per))
    found_event = list(dict.fromkeys(found_event))
    return found_per, found_event

def read_branch_sample(tree, branch, max_entries=10):
    try:
        arr = tree[branch].array(library="np", entry_stop=max_entries)
        return arr
    except Exception as e:
        print(f"  couldn't read branch {branch}: {e}")
        return None

def is_jagged_array(arr):
    # uproot returns numpy arrays; jagged arrays come as awkward arrays
    import awkward as ak
    return ak.is_awkward(arr) and arr.ndim == 1

def extract_nu_list(tree, nevents=None):
    """
    尝试提取投射侧（projectile-side）每个核子的碰撞次数 n_i 的逻辑。
    对 p+Au (NpartA==1) 的特殊处理：若存在 event-level Ncoll scalar，则用它作为 projectile nu。
    返回:
      - nu_list: for each projectile nucleon (for pA, one entry per event) an integer nu
      - meta: dict with info what method was used
    """
    # detect candidate branches
    found_per, found_event = detect_candidate_branches(tree)
    print("Detected candidate per-nucleon branches:", found_per)
    print("Detected candidate event-level branches:", found_event)

    # Attempt 1: look for obvious jagged arrays per projectile-side, e.g. 'NcollA' as jagged
    for cand in found_per:
        # read small sample to see structure
        try:
            arr = tree[cand].array(library="ak")  # awkward
            import awkward as ak
            # if it's jagged (array of arrays), check shape
            if isinstance(arr, ak.Array):
                # we'll assume it's per-nucleon list; but need to know which side (A or B)
                print(f"Using jagged per-nucleon branch: {cand}")
                # flatten for projectile side selection later; but we need per-event lists
                # Decide which side this array corresponds to by name heuristics:
                side = "A" if "a" in cand.lower() else ("B" if "b" in cand.lower() else "A")
                nu_per_event = []
                for sub in arr:
                    # if sub is empty, skip
                    if len(sub) == 0:
                        nu_per_event.append([])
                    else:
                        # ensure integers
                        nu_per_event.append([int(x) for x in sub.tolist()])
                return nu_per_event, {"method":"jagged_branch", "branch":cand, "side":side}
        except Exception as e:
            # continue to next candidate
            continue

    # Attempt 2: p+A, if NpartA exists and equals ~1 for many events, and event-level Ncoll exists:
    # Check presence of NpartA
    available = list(tree.keys())
    lower_av = [n.lower() for n in available]
    # prefer possible scalar names 'NcollA' or 'fNcollA' etc.
    scalar_candidates = [b for b in found_event if 'ncoll' in b.lower()]
    # read npartA if possible
    npartA_name = None
    for nm in available:
        if nm.lower() in ("nparta","npart_a","fNpartA".lower()):
            npartA_name = nm
            break
    # fallback: search for 'NpartA' in any case-insensitive way
    if npartA_name is None:
        for nm in available:
            if nm.lower().startswith("nparta"):
                npartA_name = nm
                break

    # read small samples to inspect if NpartA is scalar
    if npartA_name and npartA_name in tree.keys():
        try:
            sample_npartA = tree[npartA_name].array(entry_stop=100, library="np")
            mean_npartA = float(np.mean(sample_npartA))
            print(f"Mean {npartA_name} (first 100 events) = {mean_npartA:.3f}")
        except Exception as e:
            mean_npartA = None
    else:
        mean_npartA = None

    # If looks like p+A (mean npartA ~ 1), then try to use scalar Ncoll as projectile nu
    if mean_npartA is not None and mean_npartA < 2.1:
        # find the best scalar Ncoll candidate
        chosen = None
        for cand in scalar_candidates:
            # read first 10 to make sure scalar
            try:
                arr = tree[cand].array(entry_stop=10, library="np")
                if np.asarray(arr).ndim == 1:
                    chosen = cand
                    break
            except Exception:
                continue
        if chosen:
            print(f"Using event-level scalar branch {chosen} as projectile nu (p+A fallback).")
            # read full array
            arr_full = tree[chosen].array(library="np")
            # return list of single-element lists to be consistent: one projectile nucleon per event
            nu_per_event = [[int(x)] for x in arr_full.tolist()]
            return nu_per_event, {"method":"scalar_Ncoll_as_projectile_nu", "branch":chosen}

    # Attempt 3: final fallback: if there's any scalar Ncoll, use it (but warn)
    for cand in scalar_candidates:
        try:
            arr = tree[cand].array(entry_stop=10, library="np")
            if np.asarray(arr).ndim == 1:
                print(f"Fallback: using scalar branch {cand} as event-level total Ncoll (caveat).")
                arr_full = tree[cand].array(library="np")
                nu_per_event = [[int(x)] for x in arr_full.tolist()]
                return nu_per_event, {"method":"fallback_scalar_Ncoll", "branch":cand}
        except Exception:
            continue

    raise RuntimeError("Couldn't detect per-nucleon nor reasonable event-level Ncoll branch. "
                       "Please inspect branches (script printed them) and let me know the branch names for per-nucleon counts.")

def sample_final_rapidity_for_event(n_list, yb, alpha1, alpha_rest, rng):
    """
    n_list: list of integers (collision counts) for nucleons of the projectile side in one event.
            for p+A this is typically [nu] (single-element).
    returns list of final rapidities for those projectile nucleons.
    """
    y_final_list = []
    for n in n_list:
        y = yb
        n_int = int(n)
        if n_int > 0:
            dy = rng.exponential(1.0/alpha1)
            y -= dy
            if n_int > 1:
                dy_rest = rng.exponential(1.0/alpha_rest, size=n_int-1)
                y -= dy_rest.sum()
        # if n==0, final y is yb (no collisions)
        y_final_list.append(y)
    return y_final_list

def apply_centrality_cut(tree, mask_name, low, high):
    # mask_name e.g. 'Npart' or 'B'
    arr = tree[mask_name].array(library="np")
    mask = (arr >= low) & (arr <= high)
    return mask

def main():
    print("Opening file:", FILENAME)
    tree = list_tree_branches(FILENAME, TREENAME)

    # detect & extract nu list per event
    nu_per_event, meta = extract_nu_list(tree)
    print("Extraction meta:", meta)
    n_events = len(nu_per_event)
    print(f"Number of events processed (nu_per_event length): {n_events}")

    # optional centrality cut
    if CENTRALITY_CUT:
        key, low, high = CENTRALITY_CUT
        print(f"Applying centrality cut {key} in [{low},{high}]")
        arr_key = tree[key].array(library="np")
        mask = (arr_key >= low) & (arr_key <= high)
        # filter both nu_per_event and arr_key
        nu_per_event = [nu_per_event[i] for i in range(len(nu_per_event)) if mask[i]]
        print("After cut, events:", len(nu_per_event))

    # Monte Carlo sampling of final rapidities
    rng = np.random.default_rng(seed=12345)
    y_finals = []
    for nlist in nu_per_event:
        ylist = sample_final_rapidity_for_event(nlist, YBEAM, ALPHA1, ALPHAREST, rng)
        # for pA typically single entry per event; append all
        y_finals.extend(ylist)

    y_finals = np.array(y_finals)
    print("Sampled final rapidities (summary): n=", y_finals.size,
          " mean=", np.mean(y_finals) if y_finals.size>0 else float('nan'),
          " std=", np.std(y_finals) if y_finals.size>0 else float('nan'))

    # histogram and plot
    hist, edges = np.histogram(y_finals, bins=NBINS, range=YRANGE)
    centers = 0.5*(edges[:-1]+edges[1:])

    plt.figure(figsize=(8,5))
    plt.step(centers, hist, where='mid', label=f"model alpha1={ALPHA1}, alpha_rest={ALPHAREST}")
    plt.xlabel("y")
    plt.ylabel("counts (arbitrary)")
    plt.title("Final projectile nucleon rapidity distribution (model)")
    plt.legend()
    plt.grid(True)
    outpng = os.path.splitext(os.path.basename(FILENAME))[0] + "_model_dNdy.png"
    plt.savefig(outpng, dpi=150)
    print("Saved histogram to", outpng)
    plt.show()

if __name__ == "__main__":
    main()
