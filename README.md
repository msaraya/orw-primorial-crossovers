# ORW Primorial Crossovers

Repository for the paper:

**Graded All-Orders Asymptotics for Successive-Primorial Crossovers in the Odlyzko–Rubinstein–Wolf Exponential-Integral Approximation**

Author: **Mohamed Hassan Ali Awaad Saraya**

Public preprint (Zenodo, v1): https://zenodo.org/records/22132995

## Contents

- `paper/successive_primorial_transitions_JNT_submission.pdf` — public preprint / journal-submission manuscript.
- `code/numerical_check.py` — ancillary numerical consistency check used to reproduce the numerical table in the manuscript.

## Scope

The unconditional results concern the explicitly defined Odlyzko–Rubinstein–Wolf (ORW) exponential-integral approximation. The transfer discussion for actual consecutive-prime gaps is stated separately and conditionally in the paper.

## Reproducing the numerical check

Requires Python 3 and `mpmath`:

```bash
python -m pip install mpmath
python code/numerical_check.py
```

The numerical computation is ancillary; none of the proofs depends on it.

## License

The public paper is released on Zenodo under **CC BY-NC-ND 4.0**. See the Zenodo record for the authoritative deposit metadata and license terms.
