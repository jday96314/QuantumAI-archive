"""Use the unmodified metric readers and score on the current best package."""
import importlib.util
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
reference=HERE.parent.parent/'5/integration/audit_wave3_metric.py'
spec=importlib.util.spec_from_file_location('round5_metric_audit',reference)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
audit.OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE/'best'
audit.main()
