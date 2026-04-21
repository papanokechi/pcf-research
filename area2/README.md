# Self-Adjoint Structure of PCF Recurrences

Supplementary code and data for:

**"Self-Adjoint Structure of Polynomial Continued Fraction Recurrences and its Arithmetic Consequences"**

papanokechi, submitted to *Journal of Differential Equations*, 2026.

## Key Results

- **General theorem**: For any PCF with a(n)=1 and b(n) ∈ ℤ[n], the Wallis recurrence produces the exact ODE d/dx[b(x)y'] + c(x)y = 0.
- **Apparent singularities**: Every simple root of b(x) is an apparent singularity (indicial exponents {0,0}, trivial monodromy).
- **Discriminant universality**: The CM discriminant in PCF analyses equals disc(b(x)), the discriminant of the Wallis characteristic polynomial.
- **Perron–Bessel**: PCF with b(n)=2(n+1) has limit I₀(1)/I₁(1), verified to 55 decimal digits.

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/area2_general_theory.py` | Proves and verifies the general B=A' condition symbolically (sympy) and numerically (mpmath). Covers d=1,2,3, counterexample, Sturm–Liouville weight, discriminant theorem. |
| `scripts/area2_stieltjes_test.py` | Tests and refutes the naive Stieltjes integral hypothesis. Calibration: b(n)=2n+1 gives coth(1), b(n)=2(n+1) gives I₀(1)/I₁(1). |

## Reproduce

### General theorem verification

```bash
python scripts/area2_general_theory.py
```

### Stieltjes hypothesis test

```bash
python scripts/area2_stieltjes_test.py
```

### Requirements

Python 3.10+, mpmath >= 1.3.0, sympy >= 1.12

## Key numerical certificates

| Result | Value | Digits |
|--------|-------|--------|
| I₀(1)/I₁(1) via PCF b(n)=2(n+1) | 2.24019372387... | 55 |
| I₀(2)/I₁(2) via PCF b(n)=n+1 | 1.43312742672... | 30 |
| coth(1) via PCF b(n)=2n+1 | 1.31303528719... | 50 |
| V_quad apparent singularity p₀ | 1 exactly | symbolic |
| V_quad apparent singularity q₀ | 0 exactly | symbolic |

## Connection to V_quad

The V_quad constant (see [`vquad/`](../vquad/) directory) is the degree-2 application of the general theory: b(x) = 3x² + x + 1, disc(b) = −11, both roots apparent. See also Nonlinearity submission NON-110708.
