"""Scan all .kmb and dfthcsos candidate directories to check their stats and whether they are valid full instances."""

from pathlib import Path
import json

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")

def check_all():
    print("Checking available full-candidate artifacts...")
    # Look for invsqrt
    for p in ROOT.glob("initial_experiments/**/invsqrt*.kmb"):
        print(f"invsqrt: {p} ({p.stat().st_size} bytes)")
    # Look for dqi
    for p in ROOT.glob("initial_experiments/**/dqi*.kmb"):
        print(f"dqi: {p} ({p.stat().st_size} bytes)")
    # Look for dfthcsos
    for p in ROOT.glob("initial_experiments/**/dfthcsos"):
        if p.is_dir() and (p / "submission.npz").exists():
            print(f"dfthcsos: {p}")

if __name__ == "__main__":
    check_all()
