#!/usr/bin/env python3
import uproot
import numpy as np
import matplotlib.pyplot as plt

def extract_nu_distribution(filename, treename="nt_p_Au"):
    # 打开 ROOT 文件
    with uproot.open(filename) as f:
        if treename not in f:
            raise KeyError(f"TTree '{treename}' not found in {filename}")
        tree = f[treename]
        # pA情形下, 每个事件的Ncoll即为该proton的碰撞次数ν
        nu = tree["Ncoll"].array(library="np")
    return nu

def plot_nu_distribution(nu, outname="nu_distribution.png"):
    mean, std, median = np.mean(nu), np.std(nu), np.median(nu)
    print(f"⟨ν⟩ = {mean:.3f}, σ = {std:.3f}, median = {median:.3f}")

    plt.figure(figsize=(7,5))
    bins = np.arange(0, max(nu)+2) - 0.5
    plt.hist(nu, bins=bins, histtype='stepfilled', alpha=0.7, color='C0', edgecolor='k')
    plt.xlabel("Number of collisions ν", fontsize=13)
    plt.ylabel("Counts", fontsize=13)
    plt.title(rf"$\langle \nu \rangle = {mean:.2f}$", fontsize=14)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outname, dpi=200)
    plt.show()

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Analyze ν-distribution from TGlauberMC p+A output.")
    parser.add_argument("rootfile", help="Path to TGlauberMC output ROOT file (e.g. pAu.root)")
    args = parser.parse_args()

    nu = extract_nu_distribution(args.rootfile)
    plot_nu_distribution(nu)

if __name__ == "__main__":
    main()
