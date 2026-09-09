#!/usr/bin/env python3
"""
Quick helper to visualize proton trajectories and a spherical cell boundary.

Usage:
  python tools/plot_trajectories.py --file /path/to/ProtonTrajectories_*.dat --radius 5e-6 --center 0 0 0 [--proton 0]

The input file is produced when running the simulator with `--record-positions`.
The file format now includes an event kind column:
  kind 0 = regular sample, 1 = hit point, 2 = post-reflection position.
"""

import argparse
import math
from collections import defaultdict
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def read_trajectories(path: str, proton_filter: int | None) -> Tuple[Dict[int, List[Tuple[float, float, float, float, int]]], List[int]]:
    trajectories: Dict[int, List[Tuple[float, float, float, float, int]]] = defaultdict(list)
    proton_ids: set[int] = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 6:
                continue
            t, pid, kind, x, y, z = parts
            pid_int = int(pid)
            if proton_filter is not None and pid_int != proton_filter:
                continue
            proton_ids.add(pid_int)
            trajectories[pid_int].append((float(t), float(x), float(y), float(z), int(kind)))
    # sort by time to ensure clean plotting
    for pid in trajectories:
        trajectories[pid].sort(key=lambda p: p[0])
    return trajectories, sorted(proton_ids)


def plot_sphere(ax, center: Tuple[float, float, float], radius: float):
    u, v = np.mgrid[0 : 2 * np.pi : 40j, 0 : np.pi : 20j]
    x = center[0] + radius * np.cos(u) * np.sin(v)
    y = center[1] + radius * np.sin(u) * np.sin(v)
    z = center[2] + radius * np.cos(v)
    ax.plot_wireframe(x, y, z, color="gray", linewidth=0.5, alpha=0.5)


def main():
    parser = argparse.ArgumentParser(description="Plot proton trajectories with spherical cell boundary")
    parser.add_argument("--file", required=True, help="Path to ProtonTrajectories_*.dat output")
    parser.add_argument("--radius", type=float, required=True, help="Cell radius (same units as file, e.g., meters)")
    parser.add_argument("--center", type=float, nargs=3, default=(0.0, 0.0, 0.0), help="Cell center (x y z)")
    parser.add_argument("--proton", type=int, default=None, help="Plot only a single proton id (optional)")
    args = parser.parse_args()

    trajectories, proton_ids = read_trajectories(args.file, args.proton)
    if not trajectories:
        raise SystemExit("No trajectories found to plot.")

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    plot_sphere(ax, tuple(args.center), args.radius)

    for pid in proton_ids:
        pts = np.array([(p[1], p[2], p[3]) for p in trajectories[pid]])
        kinds = np.array([p[4] for p in trajectories[pid]])
        ax.plot3D(pts[:, 0], pts[:, 1], pts[:, 2], label=f"proton {pid}", linewidth=2)

        # scatter all points with color by kind
        cmap = {0: "blue", 1: "orange", 2: "red"}
        labels = {0: "sample", 1: "hit", 2: "post-reflection"}
        for kind in (0, 1, 2):
            mask = kinds == kind
            if np.any(mask):
                ax.scatter(pts[mask, 0], pts[mask, 1], pts[mask, 2], s=12, color=cmap[kind], alpha=0.8, label=f"{labels[kind]} (p{pid})")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
