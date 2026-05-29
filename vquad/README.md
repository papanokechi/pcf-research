# V_quad Resurgence Analysis

Scripts and data for:
"A non-classical Painlevé V transcendent from a quadratic polynomial continued fraction: surface classification and resurgent Stokes data"
Nonlinearity, submission NON-110708 (2026). [Surface label corrected PIII(D6) -> PV / Sakai D5^(1), W(A3^(1)) per SIARC cascade VQ-N1.]

## Scripts (in order of execution)

| Script | Description | SIARC Iteration |
|--------|-------------|-----------------|
| `scripts/t2_iter17_stokes.py` | Initial Stokes multiplier computation via Borel summation | T2 iter 17 |
| `scripts/t2_iter18_painleve.py` | Painlevé V (Sakai D5^(1)) connection and Borel plane analysis | T2 iter 18 |
| `scripts/t2_iter19_resurgence.py` | Resurgence structure: Borel singularity at ξ₀ = 2/√3 | T2 iter 19 |
| `scripts/t2_iter20_stokes_constant.py` | Stokes constant S extraction (first version) | T2 iter 20 |
| `scripts/t2_iter20_stokes_constant_v2.py` | Stokes constant S extraction (improved, 8-digit) | T2 iter 20 |
| `scripts/t2_iter21_hyperasymptotic.py` | Hyperasymptotic expansion and alien derivative | T2 iter 21 |
| `scripts/t2_iter22_s_precision.py` | High-precision Stokes constant: S ≈ 0.43770528 | T2 iter 22 |
| `scripts/t2_iter23_jimbo.py` | Jimbo connection-formula search (PIII retained as negative evidence; surface PV / D5^(1)) | T2 iter 23 |
| `scripts/jimbo_final.py` | Final Jimbo connection formula implementation | T2 iter 23 |
| `scripts/t2_iter24_sigma_conn.py` | Connection parameter σ_conn ≈ 0.060877 | T2 iter 24 |
| `scripts/verify_frobenius_apparent.py` | Frobenius verification: apparent singularities, indicial exponents {0,0} | T2 iter 17 |

## Key results

- V_quad = 1.19737399068835760244...
- PV (Sakai D5^(1), symmetry W(A3^(1)) ≅ affine S₄): α = 1/6, β = γ = 0, δ = −1/2; generic non-classical, does not reduce to PIII (δ ≠ 0)
- Borel singularity: ξ₀ = 2/√3
- Branch exponent: −1/(3√3) (8 digits)
- Stokes constant: S ≈ 0.43770528 (8 digits)
- Connection parameter: σ_conn ≈ 0.060877

## Reproduce

```bash
python scripts/verify_frobenius_apparent.py
python scripts/t2_iter22_s_precision.py
python scripts/t2_iter24_sigma_conn.py
```
