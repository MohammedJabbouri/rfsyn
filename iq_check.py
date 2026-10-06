# PYTHON

from pathlib import Path
import math
import struct
import sys

folder = Path("output_phase1_check")
paths = [
    folder / "example_000000.iq",
    folder / "example_000001.iq",
]

expected_samples = 200
contents = []
failed = False

for path in paths:
    try:
        data = path.read_bytes()
    except OSError as error:
        print(f"FAIL: {path}: {error}")
        failed = True
        continue

    if len(data) != expected_samples * 8:
        print(f"FAIL: {path}: expected 1600 bytes, got {len(data)}")
        failed = True
        continue

    samples = list(struct.iter_unpack("<ff", data))
    contents.append(data)

    finite = all(math.isfinite(i) and math.isfinite(q) for i, q in samples)

    if not finite:
        print(f"FAIL: {path}: contains NaN or infinity")
        failed = True
        continue

    peak = max(math.hypot(i, q) for i, q in samples)
    nonzero = sum(i != 0.0 or q != 0.0 for i, q in samples)
    varying = any(sample != samples[0] for sample in samples[1:])
    zero_imaginary = all(q == 0.0 for _, q in samples)

    print(f"\n{path}")
    print(f"  Samples: {len(samples)}")
    print(f"  All finite: {finite}")
    print(f"  Nonzero samples: {nonzero}")
    print(f"  Peak magnitude: {peak:.9e}")
    print(f"  Nonconstant: {varying}")
    print(f"  Imaginary components zero: {zero_imaginary}")

    if peak <= 1e-9 or not varying or not zero_imaginary:
        print("  FAIL: expected a varying real-valued FDTD signal with AWGN disabled")
        failed = True

if len(contents) == 2:
    identical = contents[0] == contents[1]
    print(f"\nFiles identical: {identical}")
    if not identical:
        failed = True
else:
    failed = True

print("\nPASS" if not failed else "\nFAIL")
sys.exit(1 if failed else 0)