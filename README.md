i# ORW Primorial Crossovers

[![Zenodo](https://img.shields.io/badge/Zenodo-record%2022132995-blue)](https://zenodo.org/records/22132995)
[![Paper license: CC BY-NC-ND 4.0](https://img.shields.io/badge/Paper%20license-CC%20BY--NC--ND%204.0-lightgrey)](https://creativecommons.org/licenses/by-nc-nd/4.0/)
[![Code license: MIT](https://img.shields.io/badge/Code%20license-MIT-green)](code/LICENSE)

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

Licensing is split by material:

- The **paper and repository documentation** are licensed under **CC BY-NC-ND 4.0**.
- The **ancillary numerical-check code** under `code/` is licensed under the **MIT License**.

See `LICENSE.md` and `code/LICENSE` for details. The Zenodo record remains the authoritative public deposit for the paper.
