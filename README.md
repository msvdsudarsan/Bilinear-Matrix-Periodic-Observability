# Bilinear Matrix Periodic Systems — Observability (MATLAB Codes)

**Reproducible Research | MATLAB Implementation | Numerical Verification**

This repository contains MATLAB scripts used to reproduce the numerical experiments reported in the research paper

“Equivalence of Kalman and Hewer Observability for Generalised Bilinear Matrix Periodic Systems.”

---

## Authors

Sri Venkata Durga Sudarsan Madhyannapu
Sravanam Pradheep Kumar

---

## Overview

This repository provides a reproducible MATLAB implementation of the numerical experiments used to verify theoretical results on **observability of generalised bilinear matrix periodic systems**.

The MATLAB scripts implement the lifted vector representation of the matrix system and compute the **periodic observability Gramian**. The numerical experiments confirm the theoretical equivalence between **Kalman observability** and **Hewer observability** for this class of periodic systems.

All numerical results reported in the manuscript have been independently verified using an additional RK45-based implementation. 

---

## Mathematical Model

The matrix differential system considered in the paper is

Ẋ(t) = A(t)X(t) + X(t)B(t) + F(t)X(t)G(t) + K(t)U(t)

with output equation

Y(t) = C(t)X(t)

Using the vectorisation identity

Vec(FXG) = (Gᵀ ⊗ F) Vec(X)

the matrix system is transformed into a lifted vector system of dimension (n^2).

Observability analysis is then performed using the **periodic observability Gramian** and the spectral properties of the **monodromy matrix**.

---

## Numerical Example

The numerical verification uses a representative system with

n = 3
lifted dimension = 9
period (T = 2\pi)

---

## Monodromy Matrix Property

The minimum singular value of the monodromy matrix is

σ_min(M) = **2.2941 × 10⁻³**

This quantity is used in the spectral lower bound for the observability Gramian.

---

## Observability Gramian Statistics

For the periodic observability Gramian (O_i):

**Period i = 1**

λ_min(O₁) = 1.2807
λ_max(O₁) = 8.413 × 10⁶
κ(O₁) = 6.569 × 10⁶

**Period i = 2**

λ_min(O₂) = 1.6946
λ_max(O₂) = 1.461 × 10⁹
κ(O₂) = 8.621 × 10⁸

**Period i = 3**

λ_min(O₃) = 2.0129
λ_max(O₃) = 5.588 × 10¹⁰
κ(O₃) = 2.776 × 10¹⁰

All computed Gramians are **positive definite**, confirming observability of the lifted system. 

---

## Spectral Lower Bound (Key Result)

The paper establishes the bound

λ_min(O_i) ≥ σ_min(M)^{2(i−1)} λ_min(O₁)

Using

σ_min(M) = 2.2941 × 10⁻³

the theoretical bounds are

i = 2

lower bound = **6.740 × 10⁻⁶**
actual value λ_min(O₂) = **1.6946**

i = 3

lower bound = **3.547 × 10⁻¹¹**
actual value λ_min(O₃) = **2.0129**

The observed values are significantly larger than the theoretical bounds, confirming the validity of the result. 

---

## Additional Numerical Verification

Further numerical experiments include:

Recursion residual verification

‖O₂ − (O₁ + MᵀO₁M)‖₍F₎ = **0**

H-observability condition

minimum integral value = **17.8683 > 0**

Minimum detection energy

J*₍obs₎ = **0.3490**

These results confirm the equivalence between **Kalman observability** and **Hewer observability**. 

---

## Scaling Study

Condition numbers of the observability Gramian across different system dimensions:

n = 2

κ(O₁) ≈ **3.794 × 10⁹ ± 1.194 × 10¹⁰**

n = 3

κ(O₁) ≈ **1.243 × 10¹⁴ ± 3.874 × 10¹⁴**

Large condition numbers arise due to unstable monodromy directions, which is typical for periodic systems with strong spectral growth. 

---

## Repository Structure

```
Bilinear-Matrix-Periodic-Observability
│
├── Supp_Paper2_Observability_v3.m
│   Main MATLAB script performing the numerical verification
│
├── Generate_Figures_Paper2.m
│   Script used to generate the figures appearing in the paper
│
├── Fig1_Obs_EigModuli.pdf
├── Fig2_Obs_GramianGrowth.pdf
├── Fig3_Obs_SpectralBound.pdf
├── Fig4_Obs_ScalingStudy.pdf
│
├── LICENSE
└── README.md
```

---

## MATLAB Environment

The numerical experiments were executed using

MATLAB R2024b

Numerical integration uses

ode45

with solver tolerances

RelTol = 1e-8
AbsTol = 1e-10

---

## Reproducing the Numerical Results

Run

Supp_Paper2_Observability_v3.m

to reproduce the numerical verification.

Run

Generate_Figures_Paper2.m

to regenerate the figures used in the manuscript.

---

## Repository Link

https://github.com/msvdsudarsan/Bilinear-Matrix-Periodic-Observability

---

## License

This repository is distributed under the MIT License.
