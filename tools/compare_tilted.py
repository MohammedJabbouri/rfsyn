import argparse
import csv
import math
import subprocess
from pathlib import Path


def read_samples(filename, angle, shift, polarization, scale):
    with filename.open(newline="") as stream:
        rows = sorted(
            csv.DictReader(stream),
            key=lambda row: int(row["probe_id"]),
        )

    probe_ids = [int(row["probe_id"]) for row in rows]

    if len(rows) != 4 or probe_ids != list(range(4)):
        raise RuntimeError(
            f"Incomplete or duplicate probe data: {filename}"
        )

    for row in rows:
        if int(row["measurement_converged"]) != 1:
            raise RuntimeError(
                f"Measurement did not settle: {filename}"
            )

        metadata_matches = (
            int(row["scale"]) == scale
            and row["polarization"] == polarization
            and math.isclose(float(row["angle_deg"]), angle,
                abs_tol=1e-12,
            )
            and math.isclose(
                float(row["shift_mm"]),
                shift,
                abs_tol=1e-12,
            )
        )

        if not metadata_matches:
            raise RuntimeError(
                f"Mismatched geometry metadata: {filename}"
            )

        for key in (
            "incident_re",
            "incident_im",
            "scattered_re",
            "scattered_im",
            "window_drift",
        ):
            if not math.isfinite(float(row[key])):
                raise RuntimeError(
                    f"Nonfinite {key}: {filename}"
                )

    return rows


def scattered_vector(rows):
    return [
        complex(
            float(row["scattered_re"]),
            float(row["scattered_im"]),
        )
        for row in rows
    ]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare a fixed finite tilted slab across "
            "uniform-grid resolutions."
        )
    )

    parser.add_argument(
        "executable",
        type=Path,
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("accuracy-results/tilted"),
    )

    parser.add_argument(
        "--reference-scale",
        type=int,
        choices=(4, 8),
        default=4,
    )

    parser.add_argument(
        "--full",
        action="store_true",
        help=(
            "Run four angles + four physical offsets, "
            "and both polarizations. Expensive."
        ),
    )

    args = parser.parse_args()

    executable = args.executable.resolve()

    if not executable.is_file():
        parser.error(
            f"Executable not found: {executable}"
        )

    args.output.mkdir(
        parents=True,
        exist_ok=True,
    )

    angles = (
        (0, 15, 30, 45)
        if args.full
        else (30,)
    )

    shifts = (
        (0.0, 0.25, 0.5, 0.75)
        if args.full
        else (0.5,)
    )

    polarizations = (
        ("TE", "TM")
        if args.full
        else ("TE",)
    )

    scales = [
        scale
        for scale in (1, 2, 4, 8)
        if scale <= args.reference_scale
    ]

    comparisons = []

    for angle in angles:
        for shift in shifts:
            for polarization in polarizations:
                runs = {}

                for scale in scales:
                    case_name = (
                        f"a{angle}_shift{shift:g}_"
                        f"{polarization}_s{scale}"
                    )

                    case_directory = (
                        args.output / case_name
                    ).resolve()

                    case_directory.mkdir(
                        parents=True,
                        exist_ok=True,
                    )

                    command = [
                        str(executable),
                        "--angle",
                        str(angle),
                        "--shift-mm",
                        str(shift),
                        "--pol",
                        polarization,
                        "--scale",
                        str(scale),
                        "--output",
                        "samples.csv",
                    ]

                    print(
                        f"Running {case_name}",
                        flush=True,
                    )

                    with (
                        case_directory / "run.log"
                    ).open("w") as log:
                        subprocess.run(
                            command,
                            cwd=case_directory,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                            check=True,
                        )

                    runs[scale] = read_samples(
                        case_directory / "samples.csv",
                        angle,
                        shift,
                        polarization,
                        scale,
                    )

                reference = runs[
                    args.reference_scale
                ]

                reference_scattered = (
                    scattered_vector(reference)
                )

                incident_norm = math.sqrt(
                    sum(
                        float(row["incident_re"]) ** 2
                        + float(row["incident_im"]) ** 2
                        for row in reference
                    )
                )

                scattered_norm = math.sqrt(
                    sum(
                        abs(value) ** 2
                        for value in reference_scattered
                    )
                )

                if incident_norm <= 1e-12:
                    raise RuntimeError(
                        "Reference illumination is too small"
                    )

                for scale, rows in runs.items():
                    if scale == args.reference_scale:
                        continue

                    for row, reference_row in zip(
                        rows,
                        reference,
                    ):
                        for key in (
                            "frequency_hz",
                            "x_m",
                            "y_m",
                        ):
                            if not math.isclose(
                                float(row[key]),
                                float(reference_row[key]),
                                rel_tol=1e-10,
                                abs_tol=1e-12,
                            ):
                                raise RuntimeError(
                                    "Incompatible "
                                    f"{key} across resolutions"
                                )

                    scattered = scattered_vector(rows)

                    difference = math.sqrt(
                        sum(
                            abs(value - target) ** 2
                            for value, target in zip(
                                scattered,
                                reference_scattered,
                            )
                        )
                    )

                    incident_normalized_difference = (
                        difference / incident_norm
                    )

                    relative_scattered_difference = (
                        difference / scattered_norm
                        if scattered_norm > 1e-12
                        else "undefined"
                    )

                    row = {
                        "angle_deg": angle,
                        "shift_mm": shift,
                        "polarization": polarization,
                        "scale": scale,
                        "reference_scale": args.reference_scale,
                        "scattered_difference_over_incident": (
                            incident_normalized_difference
                        ),
                        "scattered_relative_difference": (
                            relative_scattered_difference
                        ),
                        "window_drift": max(
                            float(item["window_drift"])
                            for item in rows
                        ),
                    }

                    comparisons.append(row)

                    print(
                        f"{polarization}, angle={angle}, "
                        f"shift={shift:g} mm, "
                        f"scale={scale} vs "
                        f"{args.reference_scale}: "
                        "scattered difference / incident = "
                        f"{incident_normalized_difference:.6g}",
                        flush=True,
                    )

    comparison_file = (
        args.output / "tilted_comparison.csv"
    )

    with comparison_file.open(
        "w",
        newline="",
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(comparisons[0]),
        )

        writer.writeheader()
        writer.writerows(comparisons)

    print(
        "MEASUREMENT COMPARISON ONLY: the finite "
        "fine-grid reference is not exact ground truth. "
        "No staircasing fix is claimed."
    )


if __name__ == "__main__":
    main()