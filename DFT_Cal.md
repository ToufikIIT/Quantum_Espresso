# 🧪 DFT Calculations for Simple Solids

---

## 📚 Table of Contents

1. [Introduction](#introduction)
2. [Periodic Structures, Supercells, and Lattice Parameters](#periodic-structures-supercells-and-lattice-parameters)
3. [Face-Centered Cubic Materials](#face-centered-cubic-materials)
4. [Hexagonal Close-Packed Materials](#hexagonal-close-packed-materials)
---

## 🔹 Introduction

In DFT studies of metals, a fundamental question is:
> **"What is the crystal structure of the lowest energy?"**

DFT helps answer this by calculating total energies for different configurations.

---

## 🔹 Periodic Structures, Supercells, and Lattice Parameters

A **simple cubic lattice** consists of atoms at:
\[
\mathbf{r} = (n_1 a, n_2 a, n_3 a)
\]

We define:
- **Unit cell**: basic repeating unit
- **Supercell**: expanded cell with periodic images
- **Lattice parameter \( a \)**: edge length

### 📈 DFT Energy vs Lattice Constant (Simple Cubic)

![Simple Cubic Energy Graph](./svgs/simple_cubic_energy.svg)

#### Taylor Expansion:
\[
E_{tot}(a) \approx E_{tot}(a_0) + \alpha(a - a_0) + \beta(a - a_0)^2
\]

#### Quadratic Fit:
\[
E_{tot}(a) = E_0 + \beta (a - a_0)^2
\]

#### Birch-Murnaghan Equation (2.3):
\[
E_{tot}(a) = E_0 + \frac{9V_0 B_0}{16} \left\{ \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^3 B_0' + \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^2 \left[ 6 - 4 \left( \frac{a_0}{a} \right)^2 \right] \right\}
\]

---

## 🔹 Face-Centered Cubic Materials (FCC)

FCC is a common structure for metals like Cu. The FCC cell has atoms at:
\[
(0,0,0), \quad (a/2,a/2,0), \quad (a/2,0,a/2), \quad (0,a/2,a/2)
\]

**Primitive cell vectors**:
\[
\mathbf{a}_1 = a\left(\frac{1}{2}, 0, \frac{1}{2}\right), \quad 
\mathbf{a}_2 = a\left(0, \frac{1}{2}, \frac{1}{2}\right), \quad 
\mathbf{a}_3 = a\left(\frac{1}{2}, \frac{1}{2}, 0\right)
\]

- Nearest neighbor distance = \( a/\sqrt{2} \)
- Cell vectors are non-orthogonal.

### 📐 FCC Lattice Geometry

![FCC Structure](./svgs/fcc_structure.svg)

### 📈 Energy vs Lattice Constant (FCC Cu)

![FCC Energy Graph](./svgs/fcc_energy.svg)

Minimum energy at:
- \( a = 3.64 \, \text{Å} \)
- \( B_0 = 142 \, \text{GPa} \)

---

## 🔹 Hexagonal Close-Packed (HCP) Materials

HCP is another densely packed structure.

**Lattice vectors:**
\[
\mathbf{a}_1 = (a, 0, 0), \quad 
\mathbf{a}_2 = \left(\frac{a}{2}, \frac{a\sqrt{3}}{2}, 0\right), \quad 
\mathbf{a}_3 = (0, 0, c)
\]

Atoms placed at:
\[
(0,0,0), \quad \left(\frac{2a}{3}, \frac{2a}{\sqrt{3}}, \frac{c}{2}\right)
\]

### 📐 HCP Lattice Geometry

![HCP Structure](./svgs/hcp_structure.svg)

### 📈 Energy vs c/a Ratio (HCP)

![HCP Energy vs c/a](./svgs/hcp_energy.svg)

> Minimum occurs near \( c/a = 1.60 \), close to ideal \( 1.633 \), but not lower than FCC → Cu is FCC in reality.

---

## 🔹 Crystal Structure Prediction

To predict structure:
- Compare DFT energies of multiple lattice types
- Use high-precision convergence
- Adjust all variables (lattice constants, atom positions)

> **HCP needs optimization of 2 parameters \( a, c \), while cubic needs 1.**

---

## ✅ Summary

- DFT accurately predicts stable lattice structures
- Simple cubic: rarely real but educational
- FCC: energetically favored for Cu
- HCP: used for comparison, same density as FCC
- Birch-Murnaghan fit improves accuracy over wide \( a \) range

---

## 📎 References

- *Density Functional Theory: A Practical Introduction* – Sholl & Steckel
- Quantum ESPRESSO Docs: https://www.quantum-espresso.org/

---

