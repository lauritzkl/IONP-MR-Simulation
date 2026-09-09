#!/usr/bin/env python3
import argparse
import struct
import numpy as np
import matplotlib.pyplot as plt

def load_bgrid(path):
    with open(path, "rb") as f:
        header = f.read(12)
        if len(header) != 12:
            raise ValueError("File too small to contain header")
        nx, ny, nz = struct.unpack("<III", header)  # little-endian uint32
        data = np.fromfile(f, dtype=np.float64)
    expected = nx * ny * nz
    if data.size != expected:
        raise ValueError(f"Expected {expected} values, got {data.size}")
    return data.reshape((nz, nx, ny)), (nx, ny, nz)  # z-major as written

def middle_indices(nz, count=10):
    if nz <= count:
        return list(range(nz))
    start = (nz - count) // 2
    return list(range(start, start + count))

def plot_slices(bgrid, output=None, vmax=None, count=10):
    nz, nx, ny = bgrid.shape
    sel = middle_indices(nz, count)
    n = len(sel)
    ncols = 5
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows), constrained_layout=True)
    vlim = 0.1
    for idx, ax in enumerate(axes.flat):
        if idx < n:
            z = sel[idx]
            im = ax.imshow(bgrid[z], origin="lower", cmap="coolwarm", vmin=-vlim, vmax=vlim)
            ax.set_title(f"z slice {z}")
        else:
            ax.axis("off")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="B field")
    if output:
        fig.savefig(output, dpi=150)
    else:
        plt.show()

def plot_compare(grid_a, grid_b, output=None, count=10):
    if grid_a.shape != grid_b.shape:
        raise ValueError(f"Grid shapes differ: {grid_a.shape} vs {grid_b.shape}")
    diff = grid_a - grid_b
    nz, nx, ny = grid_a.shape
    vlim = max(np.max(np.abs(grid_a)), np.max(np.abs(grid_b)), 1e-12)
    dlim = np.max(np.abs(diff))

    sel = middle_indices(nz, count)
    n = len(sel)
    ncols = 3
    nrows = n
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), constrained_layout=True)
    for row, z in enumerate(sel):
        axes[row, 0].imshow(grid_a[z], origin="lower", cmap="coolwarm", vmin=-0.1, vmax=0.1)
        axes[row, 0].set_title(f"A z={z}")
        axes[row, 1].imshow(grid_b[z], origin="lower", cmap="coolwarm", vmin=-0.1, vmax=0.1)
        axes[row, 1].set_title(f"B z={z}")
        axes[row, 2].imshow(diff[z], origin="lower", cmap="bwr", vmin=-dlim, vmax=dlim)
        axes[row, 2].set_title(f"Diff z={z}")
    fig.colorbar(plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(vmin=-0.1, vmax=0.1)),
                 ax=axes[:, :2].ravel().tolist(), shrink=0.5, label="B field")
    fig.colorbar(plt.cm.ScalarMappable(cmap="bwr", norm=plt.Normalize(vmin=-dlim, vmax=dlim)),
                 ax=axes[:, 2].ravel().tolist(), shrink=0.5, label="B diff")
    if output:
        fig.savefig(output, dpi=150)
    else:
        plt.show()
    return diff

def main():
    ap = argparse.ArgumentParser(description="Plot or compare B-field grids (z-major) from binary dumps.")
    ap.add_argument("path", help="Path to bgrid binary (from --save-bgrid)")
    ap.add_argument("path2", nargs="?", help="Optional second grid to compare against the first")
    ap.add_argument("--output", help="Save plot instead of showing (e.g., out.png)")
    ap.add_argument("--vmax", type=float, help="Fixed color scale max (symmetric about zero) for single-grid mode")
    args = ap.parse_args()

    bgrid, dims = load_bgrid(args.path)
    print(f"Loaded A with dimensions (nx, ny, nz)={dims}")

    if args.path2:
        bgrid_b, dims_b = load_bgrid(args.path2)
        print(f"Loaded B with dimensions (nx, ny, nz)={dims_b}")
        if dims != dims_b:
            raise ValueError(f"Grid dimensions differ: {dims} vs {dims_b}")
        diff = plot_compare(bgrid, bgrid_b, output=args.output, count=10)
        print(f"Diff stats: max={np.max(diff):.6g}, max_abs={np.max(np.abs(diff)):.6g}, "
              f"mean={np.mean(diff):.6g}, mean_abs={np.mean(np.abs(diff)):.6g}, "
              f"rmse={np.sqrt(np.mean(diff**2)):.6g}")
    else:
        plot_slices(bgrid, output=args.output, vmax=args.vmax, count=10)

if __name__ == "__main__":
    main()
