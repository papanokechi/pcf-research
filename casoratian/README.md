# Casoratian Identities — Companion to the April-10 PCF Paper

This subdirectory accompanies the manuscript

> **Papanokechi (2026).**
> *Polynomial Continued Fractions: a Proved Logarithmic Ladder,
> a 4/π Casoratian Identity, and 482 Irrational Constants.*
> Zenodo. https://doi.org/10.5281/zenodo.19491767

The full reproducibility scaffold — verification scripts (`pslq_search.py`,
`casoratian_verifier.py`), the 482-entry quadratic catalogue
(`new_irrational_constants.{json,csv}`), and dependency pin —
is maintained in a dedicated companion repository:

→ **https://github.com/papanokechi/pcf-casoratian-identities**

The split is intentional. `pcf-casoratian-identities` is a lightweight,
endorsement-ready package focused exclusively on the April-10 paper's
proofs and certificates. `pcf-research/` houses the broader SIARC stack
(area2, channel, pcf2, vquad). New users seeking only the
Casoratian-paper verification scripts should clone the dedicated repo.

## What you find here

This directory is intentionally minimal — a single README pointer.
No scripts or data are duplicated. Cloning the dedicated repo above
gives you everything referenced in the paper's Data Availability section.

## Citation

If you use these scripts or constants, please cite the Zenodo record:

```bibtex
@misc{papanokechi2026casoratian,
  author       = {Papanokechi},
  title        = {Polynomial Continued Fractions: a Proved Logarithmic Ladder,
                  a 4/$\pi$ Casoratian Identity, and 482 Irrational Constants},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19491767},
  url          = {https://doi.org/10.5281/zenodo.19491767}
}
```

(Concept DOI `10.5281/zenodo.19491767` resolves to the latest version;
v1 deposited 2026-04-10 has version-DOI `10.5281/zenodo.19491768`.)
