"""Inspect files in initial_experiments/3/submission.zip"""

import zipfile
from pathlib import Path

zpath = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/submission.zip")
with zipfile.ZipFile(zpath, "r") as z:
    for info in z.infolist():
        print(f"{info.filename}: {info.file_size} bytes")
