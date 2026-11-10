"""Analyze exact per-iteration active positions and halting behavior across all RS trajectories."""

from __future__ import annotations

import kickmix as km
import numpy as np
from verify_submission_dqi import optimal_dqi_weights


def main() -> None:
    w_probs = optimal_dqi_weights(4095, 35, 12, 2016) ** 2
    for k in range(36):
        if w_probs[k] > 1e-12 or k < 5 or k > 30:
            print(f"weight k={k:2d}: prob = {w_probs[k]:.3e}")


if __name__ == "__main__":
    main()
