#!/usr/bin/env python3
"""Run synthetic examples and validate their public output contract (stdlib only)."""

import argparse
import math
from pathlib import Path
import struct
import subprocess


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def rows(path: Path, columns: int) -> list[list[float]]:
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        require(len(values) == columns and all(map(math.isfinite, values)),
                f"Invalid numeric record in {path.name}")
        result.append(values)
    require(bool(result), f"Empty output: {path.name}")
    return result


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)


def check_signal(path: Path) -> list[list[float]]:
    signal = rows(path, 2)
    require(close(signal[0][0], 1) and close(signal[0][1], 0),
            f"Unexpected initial signal in {path.name}")
    require(all(row[0] >= 0 for row in signal), f"Negative magnitude in {path.name}")
    require(all(a[1] < b[1] for a, b in zip(signal, signal[1:])),
            f"Non-increasing signal times in {path.name}")
    return signal


def check_run(name: str, directory: Path, expected_positions: list[list[float]]) -> dict[str, bytes]:
    data = directory / "data"
    log = "\n".join(p.read_text(encoding="utf-8", errors="replace")
                    for p in (directory / "logs").glob("*.txt"))
    require("finished computation" in log, f"{name}: missing completion log")
    require("Simulation failed for" not in log and "Fatal error:" not in log,
            f"{name}: simulation failure in log")
    positions = list(data.glob("*_IONPs_*.dat"))
    require(len(positions) == 1, f"{name}: expected one IONP output")
    actual_positions = rows(positions[0], 3)
    require(actual_positions == expected_positions, f"{name}: geometry differs from input")
    sequence = {"fid": "FID", "se": "SE", "mse": "MSE", "cell_importance": "FID"}[name]
    signal_files = list(data.glob(f"*_signal_{sequence}_*.dat"))
    require(len(signal_files) == (2 if name == "mse" else 1),
            f"{name}: unexpected signal-file count")
    signals = sorted((check_signal(p) for p in signal_files), key=len)
    full = signals[-1]
    require(len(full) == 11 and close(full[-1][1], 1e-4),
            f"{name}: unexpected full signal length/end time")
    if name == "mse":
        sampled = signals[0]
        require(len(sampled) == 5 and all(close(row[1], k * 2e-5)
                                        for k, row in enumerate(sampled)),
                "mse: unexpected sampled times")
        require(all(all(close(a, b) for a, b in zip(row, full[2 * k]))
                    for k, row in enumerate(sampled)),
                "mse: sampled signal differs from full signal")
    if name == "cell_importance":
        require("IONP populations (intracellular / extracellular): 1 / 2" in log,
                "cell_importance: unexpected IONP populations")
        contributions = []
        for proposal in ("uniform", "cell_inside", "cell_near"):
            for kind in ("mean", "mis_contribution"):
                matches = list(data.glob(f"*_proposal_{proposal}_{kind}_*.dat"))
                require(len(matches) == 1, f"Missing {proposal} {kind} output")
                curve = rows(matches[0], 3)
                require(len(curve) == len(full) and all(close(a[2], b[1])
                        for a, b in zip(curve, full)), "Proposal time grid differs")
                if kind == "mis_contribution":
                    contributions.append(curve)
        for index, (magnitude, _) in enumerate(full):
            reconstructed = abs(sum(complex(curve[index][0], curve[index][1])
                                    for curve in contributions))
            require(close(reconstructed, magnitude), "MIS contributions do not reconstruct signal")
    trajectory_files = list(data.glob("*_ProtonTrajectories_*.dat"))
    require(len(trajectory_files) == 1, f"{name}: missing trajectory output")
    trajectory = rows(trajectory_files[0], 6)
    expected_count = 128 if name == "cell_importance" else 24
    require({row[1] for row in trajectory} == set(range(expected_count)),
            f"{name}: trajectory proton IDs do not match sample count")
    require(all(row[2] in (0, 1, 2) for row in trajectory), "Invalid trajectory event")
    grid = (data / "bgrid.bin").read_bytes()
    require(len(grid) == 12 + 12 ** 3 * 8, f"{name}: grid size differs")
    require(struct.unpack("=III", grid[:12]) == (12, 12, 12), "Unexpected grid dimensions")
    require(all(math.isfinite(item[0]) for item in struct.iter_unpack("=d", grid[12:])),
            "Non-finite grid value")
    return {p.name: p.read_bytes() for p in data.iterdir() if p.is_file()}


def validate(executable: Path, source: Path, work: Path) -> None:
    require(executable.is_file(), f"Executable not found: {executable}")
    require(not work.exists() or (work.is_dir() and not any(work.iterdir())),
            "Choose an absent or empty work directory; existing results are never removed")
    work.mkdir(parents=True, exist_ok=True)
    examples = source / "examples"
    positions = rows(examples / "IONP_pos/positions.txt", 3)
    for name in ("fid", "se", "mse", "cell_importance"):
        repeats = []
        for repeat in (1, 2):
            directory = work / name / str(repeat)
            directory.mkdir(parents=True)
            args = [str(executable), "--ini-dir", str(examples), "--ini-root", str(examples),
                    "--config", str(examples / f"{name}.txt"),
                    "--output-dir", str(directory / "data"),
                    "--log-dir", str(directory / "logs"), "--seed", "123",
                    "--record-positions", "--save-region-signals", "--save-bgrid",
                    str(directory / "data/bgrid.bin")]
            run = subprocess.run(args, cwd=source, capture_output=True, text=True,
                                 encoding="utf-8", errors="replace", timeout=60)
            (directory / "stdout.txt").write_text(run.stdout, encoding="utf-8")
            (directory / "stderr.txt").write_text(run.stderr, encoding="utf-8")
            require(run.returncode == 0 and "Simulation failed for" not in run.stderr
                    and "Fatal error:" not in run.stderr, f"{name}: simulation failed: {run.stderr}")
            repeats.append(check_run(name, directory, positions))
        require(repeats[0] == repeats[1], f"{name}: seeded numerical outputs differ")
        print(f"{name}: signals, geometry, grid, trajectories and seeded repeat passed", flush=True)
    print("All four publication examples passed, including MIS reconstruction.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        validate(args.executable.resolve(), args.source_dir.resolve(), args.work_dir.resolve())
    except (OSError, ValueError, subprocess.TimeoutExpired) as error:
        parser.exit(1, f"Publication example check failed: {error}\n")


if __name__ == "__main__":
    main()
