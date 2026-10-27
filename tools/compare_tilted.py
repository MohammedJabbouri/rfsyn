import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import subprocess
import sys

from datetime import datetime, timezone
from pathlib import Path


def file_hash(path):
    digest = hashlib.sha256()

    with path.open("rb") as stream:
        for block in iter(
            lambda: stream.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def atomic_csv(path, fields, rows):
    temporary = path.with_name(path.name + ".tmp")

    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
        )

        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, path)


def load_samples(path, case):
    angle, shift, polarization, scale = case

    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    rows.sort(
        key=lambda row: int(row["probe_id"])
    )

    probe_ids = [
        int(row["probe_id"])
        for row in rows
    ]

    if len(rows) != 4 or probe_ids != list(range(4)):
        raise ValueError(
            "expected exactly four distinct probes"
        )

    expected_positions = (
        (0.030, 0.032),
        (0.030, 0.040),
        (0.030, 0.048),
        (0.062, 0.040),
    )

    numeric_fields = (
        "angle_deg",
        "shift_mm",
        "dx_m",
        "frequency_hz",
        "x_m",
        "y_m",
        "incident_re",
        "incident_im",
        "scattered_re",
        "scattered_im",
        "window_drift",
        "incident_norm",
        "steps",
    )

    for row, (x, y) in zip(
        rows,
        expected_positions,
    ):
        if not all(
            math.isfinite(float(row[key]))
            for key in numeric_fields
        ):
            raise ValueError(
                "nonfinite measurement data"
            )

        checks = (
            int(row["scale"]) == scale,

            row["polarization"] == polarization,

            math.isclose(
                float(row["angle_deg"]),
                angle,
                abs_tol=1e-12,
            ),

            math.isclose(
                float(row["shift_mm"]),
                shift,
                abs_tol=1e-12,
            ),

            math.isclose(
                float(row["dx_m"]),
                0.001 / scale,
                rel_tol=1e-10,
            ),

            math.isclose(
                float(row["x_m"]),
                x,
                abs_tol=1e-12,
            ),

            math.isclose(
                float(row["y_m"]),
                y,
                abs_tol=1e-12,
            ),

            float(row["frequency_hz"]) > 0.0,

            float(row["incident_norm"]) > 0.0,

            float(row["steps"]) > 0.0,

            float(row["steps"]).is_integer(),

            int(row["measurement_converged"]) == 1,

            0.0 <= float(row["window_drift"]) <= 1e-6,
        )

        if not all(checks):
            raise ValueError(
                "geometry, resolution, or convergence "
                "check failed"
            )

        for key in (
            "frequency_hz",
            "steps",
        ):
            if not math.isclose(
                float(row[key]),
                float(rows[0][key]),
                rel_tol=1e-12,
            ):
                raise ValueError(
                    "inconsistent rows within measurement"
                )

    incident_squared = sum(
        float(row["incident_re"]) ** 2
        + float(row["incident_im"]) ** 2
        for row in rows
    )

    if incident_squared <= 1e-24:
        raise ValueError(
            "incident probe vector is too small"
        )

    return rows


def case_folder(root, case):
    angle, shift, polarization, scale = case

    return root / (
        f"a{angle:g}_shift{shift:g}_"
        f"{polarization}_s{scale}"
    )


def inspect_case(
    root,
    case,
    args,
    executable_digest,
    settings,
):
    folder = case_folder(root, case)
    samples = folder / "samples.csv"

    if not samples.is_file():
        return {
            "status": "missing",
            "note": "no samples.csv",
            "rows": None,
            "identity": "",
        }

    try:
        rows = load_samples(samples, case)

    except (
        ValueError,
        KeyError,
        TypeError,
        OSError,
        csv.Error,
    ) as error:
        return {
            "status": "invalid",
            "note": str(error),
            "rows": None,
            "identity": "",
        }

    metadata_path = folder / "run_meta.json"

    if not metadata_path.exists():
        if not args.allow_legacy_reuse:
            return {
                "status": "unverified",
                "note": (
                    "legacy data: use "
                    "--allow-legacy-reuse to accept "
                    "unknown executable provenance"
                ),
                "rows": None,
                "identity": "",
            }

        return {
            "status": "ready",
            "note": (
                "legacy: executable/settings "
                "provenance unknown"
            ),
            "rows": rows,
            "identity": "legacy_unknown",
        }

    try:
        metadata = json.loads(
            metadata_path.read_text()
        )

        if metadata["samples_sha256"] != file_hash(samples):
            raise ValueError(
                "samples hash does not match "
                "run_meta.json"
            )

        identity = metadata["executable_sha256"]

        if (
            not isinstance(identity, str)
            or len(identity) != 64
        ):
            raise ValueError(
                "invalid executable hash in metadata"
            )

        if not args.summarize_only:
            if identity != executable_digest:
                raise ValueError(
                    "executable changed; use a new "
                    "output directory, or --rerun "
                    "to replace with backups"
                )

            if metadata["settings"] != settings:
                raise ValueError(
                    "run settings changed; use a new "
                    "output directory or --rerun"
                )

    except (
        ValueError,
        KeyError,
        TypeError,
        OSError,
    ) as error:
        return {
            "status": "incompatible",
            "note": str(error),
            "rows": None,
            "identity": "",
        }

    return {
        "status": "ready",
        "note": "hash-verified saved run",
        "rows": rows,
        "identity": identity,
    }


def write_reports(
    root,
    groups,
    scales,
    reference_scale,
    records,
):
    progress_rows = []

    for case, record in records.items():
        angle, shift, polarization, scale = case

        progress_rows.append({
            "angle_deg": angle,
            "shift_mm": shift,
            "polarization": polarization,
            "scale": scale,
            "status": record["status"],
            "note": record["note"],
        })

    atomic_csv(
        root / "tilted_progress.csv",
        [
            "angle_deg",
            "shift_mm",
            "polarization",
            "scale",
            "status",
            "note",
        ],
        progress_rows,
    )

    comparisons = []

    for angle, shift, polarization in groups:
        reference = records[
            (
                angle,
                shift,
                polarization,
                reference_scale,
            )
        ]

        for scale in scales:
            if scale == reference_scale:
                continue

            run = records[
                (
                    angle,
                    shift,
                    polarization,
                    scale,
                )
            ]

            result = {
                "angle_deg": angle,
                "shift_mm": shift,
                "polarization": polarization,
                "scale": scale,
                "reference_scale": reference_scale,
                "comparison_status": "waiting_for_data",
                "provenance": "",
                "scattered_difference_over_incident": "",
                "scattered_relative_difference": "",
                "window_drift": "",
                "reference_window_drift": "",
            }

            if (
                run["status"] == "ready"
                and reference["status"] == "ready"
            ):
                known_identities = {
                    run["identity"],
                    reference["identity"],
                } - {"legacy_unknown"}

                if len(known_identities) > 1:
                    result["comparison_status"] = (
                        "blocked_mixed_executables"
                    )

                else:
                    rows = run["rows"]
                    reference_rows = reference["rows"]

                    frequency_matches = math.isclose(
                        float(rows[0]["frequency_hz"]),
                        float(
                            reference_rows[0]["frequency_hz"]
                        ),
                        rel_tol=1e-10,
                    )

                    if not frequency_matches:
                        result["comparison_status"] = (
                            "blocked_frequency_mismatch"
                        )

                    else:
                        values = [
                            complex(
                                float(row["scattered_re"]),
                                float(row["scattered_im"]),
                            )
                            for row in rows
                        ]

                        targets = [
                            complex(
                                float(row["scattered_re"]),
                                float(row["scattered_im"]),
                            )
                            for row in reference_rows
                        ]

                        incident_norm = math.sqrt(
                            sum(
                                float(row["incident_re"]) ** 2
                                + float(row["incident_im"]) ** 2
                                for row in reference_rows
                            )
                        )

                        scattered_norm = math.sqrt(
                            sum(
                                abs(value) ** 2
                                for value in targets
                            )
                        )

                        difference = math.sqrt(
                            sum(
                                abs(value - target) ** 2
                                for value, target in zip(
                                    values,
                                    targets,
                                )
                            )
                        )

                        legacy_provenance = (
                            "legacy_unknown"
                            in (
                                run["identity"],
                                reference["identity"],
                            )
                        )

                        result.update(
                            comparison_status="ready",

                            provenance=(
                                "legacy_unknown"
                                if legacy_provenance
                                else "hash_verified"
                            ),

                            scattered_difference_over_incident=(
                                difference / incident_norm
                            ),

                            scattered_relative_difference=(
                                difference / scattered_norm
                                if scattered_norm > 1e-12
                                else "undefined"
                            ),

                            window_drift=max(
                                float(row["window_drift"])
                                for row in rows
                            ),

                            reference_window_drift=max(
                                float(row["window_drift"])
                                for row in reference_rows
                            ),
                        )

            comparisons.append(result)

    comparison_fields = [
        "angle_deg",
        "shift_mm",
        "polarization",
        "scale",
        "reference_scale",
        "comparison_status",
        "provenance",
        "scattered_difference_over_incident",
        "scattered_relative_difference",
        "window_drift",
        "reference_window_drift",
    ]

    atomic_csv(
        root / "tilted_comparison.csv",
        comparison_fields,
        comparisons,
    )

    return comparisons


def run_case(
    executable,
    root,
    case,
    executable_digest,
    settings,
):
    angle, shift, polarization, scale = case

    folder = case_folder(root, case)

    folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    pending_samples = folder / "samples.pending.csv"
    pending_log = folder / "run.pending.log"

    if pending_samples.exists():
        pending_samples.unlink()

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
        pending_samples.name,

        "--settle-cycles",
        str(settings["settle_cycles"]),

        "--window-cycles",
        str(settings["window_cycles"]),

        "--max-windows",
        str(settings["max_windows"]),
    ]

    with pending_log.open("w") as log:
        child = subprocess.Popen(
            command,
            cwd=folder,
            stdout=log,
            stderr=subprocess.STDOUT,
        )

        try:
            return_code = child.wait()

        except KeyboardInterrupt:
            try:
                child.terminate()
                child.wait(timeout=5)

            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()

            except ProcessLookupError:
                pass

            raise

    if return_code != 0:
        raise RuntimeError(
            f"validator exit {return_code}; "
            f"inspect {pending_log}"
        )

    rows = load_samples(
        pending_samples,
        case,
    )

    metadata = {
        "executable_sha256": executable_digest,
        "samples_sha256": file_hash(pending_samples),
        "settings": settings,
        "command": command,
        "completed_utc": (
            datetime.now(timezone.utc).isoformat()
        ),
    }

    pending_metadata = (
        folder / "run_meta.pending.json"
    )

    pending_metadata.write_text(
        json.dumps(metadata, indent=2)
    )

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    for name in (
        "samples.csv",
        "run.log",
        "run_meta.json",
    ):
        old_file = folder / name

        if old_file.exists():
            old_file.rename(
                folder / f"{name}.bak-{timestamp}"
            )

    os.replace(
        pending_samples,
        folder / "samples.csv",
    )

    os.replace(
        pending_log,
        folder / "run.log",
    )

    os.replace(
        pending_metadata,
        folder / "run_meta.json",
    )

    return {
        "status": "ready",
        "note": "new hash-verified run",
        "rows": rows,
        "identity": executable_digest,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Resume tilted-slab runs and save "
            "partial comparisons."
        )
    )

    parser.add_argument(
        "executable",
        nargs="?",
        type=Path,
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("accuracy-results/tilted"),
    )

    parser.add_argument(
        "--full",
        action="store_true",
    )

    parser.add_argument(
        "--angles",
        nargs="+",
        type=float,
    )

    parser.add_argument(
        "--shifts",
        nargs="+",
        type=float,
        help="Physical offsets in millimeters",
    )

    parser.add_argument(
        "--pol",
        nargs="+",
        choices=("TE", "TM"),
    )

    parser.add_argument(
        "--scales",
        nargs="+",
        type=int,
        choices=(1, 2, 4, 8),
    )

    parser.add_argument(
        "--reference-scale",
        type=int,
        choices=(4, 8),
        default=4,
    )

    parser.add_argument(
        "--summarize-only",
        action="store_true",
    )

    parser.add_argument(
        "--allow-legacy-reuse",
        action="store_true",
    )

    parser.add_argument(
        "--rerun",
        action="store_true",
        help=(
            "Explicitly rerun selected cases; "
            "back up old results on success"
        ),
    )

    parser.add_argument(
        "--settle-cycles",
        type=int,
        default=24,
    )

    parser.add_argument(
        "--window-cycles",
        type=int,
        default=12,
    )

    parser.add_argument(
        "--max-windows",
        type=int,
        default=8,
    )

    args = parser.parse_args()

    angles = list(dict.fromkeys(
        args.angles
        if args.angles is not None
        else (
            [0, 15, 30, 45]
            if args.full
            else [30]
        )
    ))

    shifts = list(dict.fromkeys(
        args.shifts
        if args.shifts is not None
        else (
            [0, 0.25, 0.5, 0.75]
            if args.full
            else [0.5]
        )
    ))

    polarizations = list(dict.fromkeys(
        args.pol
        if args.pol is not None
        else (
            ["TE", "TM"]
            if args.full
            else ["TE"]
        )
    ))

    scales = sorted(set(
        args.scales
        if args.scales is not None
        else [1, 2, 4]
    ))

    if any(
        not math.isfinite(value)
        or not 0.0 <= value <= 60.0
        for value in angles
    ):
        parser.error(
            "angles must be finite and between "
            "0 and 60 degrees"
        )

    if any(
        not math.isfinite(value)
        or abs(value) > 1.0
        for value in shifts
    ):
        parser.error(
            "shifts must be finite and between "
            "-1 and 1 mm"
        )

    if any(
        scale > args.reference_scale
        for scale in scales
    ):
        parser.error(
            "scales cannot exceed the reference scale"
        )

    scales = sorted(set(
        scales + [args.reference_scale]
    ))

    if len(scales) < 2:
        parser.error(
            "select at least one scale below "
            "the reference scale"
        )

    settings_valid = (
        8 <= args.settle_cycles <= 500
        and 4 <= args.window_cycles <= 100
        and 3 <= args.max_windows <= 100
    )

    if not settings_valid:
        parser.error(
            "invalid settling/window settings "
            "for the tilted validator"
        )

    if args.summarize_only and args.rerun:
        parser.error(
            "--rerun cannot be combined "
            "with --summarize-only"
        )

    executable = (
        args.executable.resolve()
        if args.executable
        else None
    )

    if not args.summarize_only:
        if (
            executable is None
            or not executable.is_file()
        ):
            parser.error(
                "provide the validator executable, "
                "or use --summarize-only"
            )

    executable_digest = (
        file_hash(executable)
        if executable and executable.is_file()
        else None
    )

    root = args.output.resolve()

    root.mkdir(
        parents=True,
        exist_ok=True,
    )

    settings = {
        "settle_cycles": args.settle_cycles,
        "window_cycles": args.window_cycles,
        "max_windows": args.max_windows,
    }

    groups = list(itertools.product(
        angles,
        shifts,
        polarizations,
    ))

    cases = [
        (angle, shift, polarization, scale)
        for angle, shift, polarization in groups
        for scale in scales
    ]

    records = {
        case: inspect_case(
            root,
            case,
            args,
            executable_digest,
            settings,
        )
        for case in cases
    }

    if any(
        record["identity"] == "legacy_unknown"
        for record in records.values()
    ):
        print(
            "WARNING: legacy runs accepted; "
            "their executable/settings identity "
            "cannot be verified.",
            flush=True,
        )

    comparisons = write_reports(
        root,
        groups,
        scales,
        args.reference_scale,
        records,
    )

    active_case = None
    interrupted = False

    try:
        if not args.summarize_only:
            for index, case in enumerate(cases, 1):
                record = records[case]
                name = case_folder(root, case).name

                if (
                    record["status"] == "ready"
                    and not args.rerun
                ):
                    print(
                        f"[{index}/{len(cases)}] "
                        f"REUSE {name}",
                        flush=True,
                    )

                    continue

                if (
                    record["status"]
                    in ("incompatible", "unverified")
                    and not args.rerun
                ):
                    print(
                        f"[{index}/{len(cases)}] "
                        f"BLOCKED {name}: "
                        f"{record['note']}",
                        flush=True,
                    )

                    continue

                active_case = case

                print(
                    f"[{index}/{len(cases)}] "
                    f"RUN {name}",
                    flush=True,
                )

                try:
                    records[case] = run_case(
                        executable,
                        root,
                        case,
                        executable_digest,
                        settings,
                    )

                except (
                    OSError,
                    ValueError,
                    KeyError,
                    TypeError,
                    RuntimeError,
                    csv.Error,
                ) as error:
                    records[case] = {
                        "status": "failed",
                        "note": str(error),
                        "rows": None,
                        "identity": "",
                    }

                    print(
                        f"FAILED {name}: {error}",
                        flush=True,
                    )

                active_case = None

                comparisons = write_reports(
                    root,
                    groups,
                    scales,
                    args.reference_scale,
                    records,
                )

    except KeyboardInterrupt:
        interrupted = True

        if active_case is not None:
            records[active_case] = {
                "status": "interrupted",
                "note": (
                    "partial pending files "
                    "are not reusable"
                ),
                "rows": None,
                "identity": "",
            }

        print(
            "\nStopped. Completed runs and "
            "partial summaries are preserved.",
            flush=True,
        )

    finally:
        comparisons = write_reports(
            root,
            groups,
            scales,
            args.reference_scale,
            records,
        )

    ready_count = sum(
        record["status"] == "ready"
        for record in records.values()
    )

    comparison_count = sum(
        row["comparison_status"] == "ready"
        for row in comparisons
    )

    print(
        "Complete selected simulations: "
        f"{ready_count}/{len(cases)}"
    )

    print(
        "Available comparisons: "
        f"{comparison_count}/{len(comparisons)}"
    )

    print(
        "Reports: "
        f"{root / 'tilted_progress.csv'} "
        "and "
        f"{root / 'tilted_comparison.csv'}"
    )

    for row in comparisons:
        if row["comparison_status"] != "ready":
            continue

        print(
            f"{row['polarization']} "
            f"angle={row['angle_deg']:g} "
            f"shift={row['shift_mm']:g} mm "
            f"scale={row['scale']} vs "
            f"{row['reference_scale']}: "
            "difference/incident="
            f"{row['scattered_difference_over_incident']:.6g} "
            f"({row['provenance']})"
        )

    print(
        "Fine-grid comparisons are not exact "
        "physical errors or proof of a "
        "staircasing fix."
    )

    if interrupted:
        return 130

    if args.summarize_only:
        return 0

    return 0 if ready_count == len(cases) else 1


if __name__ == "__main__":
    sys.exit(main())