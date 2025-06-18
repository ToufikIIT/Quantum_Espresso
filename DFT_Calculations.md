# 📘 DFT Calculations for Simple Solids

---

## 📚 Table of Contents
1. [Introduction](#introduction)
2. [Periodic Structures, Supercells, and Lattice Parameters](#periodic-structures-supercells-and-lattice-parameters)
3. [Face-Centered Cubic Materials](#face-centered-cubic-materials)
4. [Hexagonal Close-Packed Materials](#hexagonal-close-packed-materials)
---

## 🔹 Introduction

In this chapter, we study how **Density Functional Theory (DFT)** is used to determine the **crystal structure** of solids. Though the implementation details are skipped, conceptual understanding is emphasized.

> "Our task is to define the location of all atoms in a crystal of a pure metal."

---

## 🔹 Periodic Structures, Supercells, and Lattice Parameters

- A **simple cubic metal** fills 3D space with cubes of edge length \( a \), placing a single atom at each cube corner.
- Atomic positions: \( \mathbf{r} = (n_1 a, n_2 a, n_3 a) \)

> 🟦 We define a volume that fills space: unit cell + atomic position → called a **supercell**.

### ✅ Simple Cubic Energy vs Lattice Constant

**Figure 2.1** – DFT Energy \( E_{tot} \) as function of lattice parameter \( a \):

![Simple Cubic Energy vs a](./diagrams/simple_cubic_energy.svg)

#### Equation 2.1: Taylor Expansion

\[
E_{tot}(a) \approx E_{tot}(a_0) + \alpha(a - a_0) + \beta(a - a_0)^2
\]

#### Equation 2.2: Fitted Form

\[
E_{tot}(a) = E_0 + \beta (a - a_0)^2
\]

- \( a_0 = 2.43 \, \text{Å} \), \( B_0 = 103 \, \text{GPa} \)

### Extended Fit: Birch-Murnaghan Equation (Eq 2.3)

\[
E_{tot}(a) = E_0 + \frac{9V_0 B_0}{16} \left\{ \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^3 B_0' + \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^2 \left[ 6 - 4 \left( \frac{a_0}{a} \right)^2 \right] \right\}
\]

> ✅ This allows more accurate fits over wider lattice parameter ranges.

---

## 🔹 Face-Centered Cubic (FCC) Materials

FCC is more common for metals than simple cubic.

- Atom positions:
  - \( (0,0,0), (a/2,a/2,0), (a/2,0,a/2), (0,a/2,a/2) \)
- Primitive cell vectors:
  \[
  \begin{align*}
  \mathbf{a}_1 &= a \left( \frac{1}{2}, 0, \frac{1}{2} \right) \\
  \mathbf{a}_2 &= a \left( 0, \frac{1}{2}, \frac{1}{2} \right) \\
  \mathbf{a}_3 &= a \left( \frac{1}{2}, \frac{1}{2}, 0 \right)
  \end{align*}
  \]

> 📌 Distance between nearest neighbors = \( a/\sqrt{2} \)  
> 🔺 Cell vectors are NOT orthogonal

### ✅ FCC Geometry

![FCC Structure](./diagrams/fcc_structure.svg)

### ✅ FCC Energy vs Lattice Constant

**Figure 2.3** – DFT Energy for FCC Cu vs lattice parameter:

![FCC Energy Curve](./diagrams/fcc_energy_curve.svg)

- Minimum at \( a = 3.64 \, \text{Å} \), \( B_0 = 142 \, \text{GPa} \)

> ✅ FCC Cu is more stable than simple cubic: lower energy

---

## 🔹 Hexagonal Close-Packed (HCP) Materials

> 🧩 FCC & HCP have same packing density, but HCP uses a different lattice

- Cell vectors:
  \[
  \begin{align*}
  \mathbf{a}_1 &= (a, 0, 0) \\
  \mathbf{a}_2 &= \left( \frac{a}{2}, \frac{a\sqrt{3}}{2}, 0 \right) \\
  \mathbf{a}_3 &= (0, 0, c)
  \end{align*}
  \]
- Atom positions:
  - \( (0,0,0), (2a/3, 2a/\sqrt{3}, c/2) \)

#### Equation 2.7: Fractional Coordinates

\[
\mathbf{r}_j = \sum_{j=1}^3 f_{ij} \mathbf{a}_j
\]

> 🧮 Useful to describe general atom positions in supercell via lattice vectors.

### ✅ HCP Geometry

![HCP Structure](./diagrams/hcp_structure.svg)

### ✅ HCP Total Energy vs Lattice Ratio (c/a)

**Figure 2.4** – DFT Energy vs \( a \) for different \( c/a \) ratios:

![HCP Energy Curve](./diagrams/hcp_energy_vs_ratio.svg)

> ✅ Minimum energy at \( c/a \approx 1.60 \), close to ideal \( 1.633 \). Confirms Cu is not HCP.

---

## ✅ Summary

- DFT lets us predict equilibrium structures of solids.
- Concepts like **supercells**, **fractional coordinates**, and **primitive vectors** are crucial.
- FCC structure is most stable for Cu based on energy comparison.
- Birch-Murnaghan equation gives accurate energy–lattice parameter fits.

---

## 📎 References

- *Density Functional Theory: A Practical Introduction* – Sholl & Steckel  
- Quantum ESPRESSO: [https://www.quantum-espresso.org](https://www.quantum-espresso.org)
