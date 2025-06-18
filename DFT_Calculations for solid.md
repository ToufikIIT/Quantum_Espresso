# DFT Calculations for Simple Solids

This document summarizes **Chapter 2** from *Density Functional Theory: A Practical Introduction* by Sholl and Steckel. It explains how DFT can predict the crystal structure of solids with emphasis on concepts, equations, graphs, and diagrams.

---

## 📚 Table of Contents

1. [Introduction](#introduction)
2. [Periodic Structures, Supercells, and Lattice Parameters](#periodic-structures-supercells-and-lattice-parameters)
3. [Face-Centered Cubic Materials](#face-centered-cubic-materials)
4. [Hexagonal Close-Packed Materials](#hexagonal-close-packed-materials)
5. [Crystal Structure Prediction](#crystal-structure-prediction)

---

## 🔹 Introduction

In this chapter, we study how **Density Functional Theory (DFT)** is used to determine the **crystal structure** of solids. Though the implementation details are skipped, conceptual understanding is emphasized.

> "Our task is to define the location of all atoms in a crystal of a pure metal."

---

## 🔹 Periodic Structures, Supercells, and Lattice Parameters

* A **simple cubic metal** fills 3D space with cubes of edge length $a$, placing a single atom at each cube corner.
* Atomic positions: $\mathbf{r} = (n_1 a, n_2 a, n_3 a)$

> <mark>We define a volume that fills space: unit cell + atomic position → called a **supercell**.</mark>

### ✅ Graph: Simple Cubic Energy vs Lattice Constant

**Figure 2.1** – DFT Energy $E_{tot}$ as function of lattice parameter $a$:

```svg
<svg viewBox="0 0 300 200" xmlns="http://www.w3.org/2000/svg">
  <line x1="40" y1="10" x2="40" y2="180" stroke="black"/>
  <line x1="40" y1="180" x2="290" y2="180" stroke="black"/>
  <path d="M 40 150 Q 120 100, 200 150" fill="none" stroke="blue" stroke-width="2"/>
  <text x="20" y="20">E<sub>tot</sub></text>
  <text x="240" y="195">a (Å)</text>
</svg>
```

#### Equation 2.1: Taylor Expansion

$E_{tot}(a) \approx E_{tot}(a_0) + \alpha(a - a_0) + \beta(a - a_0)^2$

#### Equation 2.2: Fitted Form

$E_{tot}(a) = E_0 + \beta (a - a_0)^2$

* $a_0 = 2.43 \, \text{Å}$, $B_0 = 103 \, \text{GPa}$

### Extended Fit: Birch-Murnaghan Equation (Eq 2.3)

$$
E_{tot}(a) = E_0 + \frac{9V_0 B_0}{16} \left\{ \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^3 B_0' + \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^2 \left[ 6 - 4 \left( \frac{a_0}{a} \right)^2 \right] \right\}
$$

> <mark>This allows more accurate fits over wider lattice parameter ranges.</mark>

---

## 🔹 Face-Centered Cubic (FCC) Materials

The FCC structure is more common than simple cubic for metals.

* Atoms placed at:

  * $(0,0,0), (a/2,a/2,0), (a/2,0,a/2), (0,a/2,a/2)$
* **Primitive cell vectors**:
  $\mathbf{a}_1 = a \left( \frac{1}{2}, 0, \frac{1}{2} \right), \quad \mathbf{a}_2 = a \left( 0, \frac{1}{2}, \frac{1}{2} \right), \quad \mathbf{a}_3 = a \left( \frac{1}{2}, \frac{1}{2}, 0 \right)$

> <mark>Distance between nearest neighbors = $a/\sqrt{2}$</mark> <mark>Cell vectors are NOT orthogonal</mark>

### ✅ FCC Geometry

```svg
<svg viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">
  <rect x="50" y="50" width="100" height="100" fill="none" stroke="black"/>
  <circle cx="50" cy="50" r="5" fill="red"/>
  <circle cx="150" cy="50" r="5" fill="red"/>
  <circle cx="50" cy="150" r="5" fill="red"/>
  <circle cx="150" cy="150" r="5" fill="red"/>
  <circle cx="100" cy="100" r="5" fill="blue"/>
</svg>
```

### ✅ Graph: FCC Energy vs Lattice Constant

**Figure 2.3** – DFT Energy for FCC Cu vs lattice parameter:

```svg
<svg viewBox="0 0 300 200" xmlns="http://www.w3.org/2000/svg">
  <line x1="40" y1="10" x2="40" y2="180" stroke="black"/>
  <line x1="40" y1="180" x2="290" y2="180" stroke="black"/>
  <path d="M 40 150 Q 150 80, 260 150" fill="none" stroke="green" stroke-width="2"/>
  <text x="20" y="20">E</text>
  <text x="250" y="195">a (Å)</text>
</svg>
```

* Minimum at $a = 3.64 \text{ Å}$, $B_0 = 142 \text{ GPa}$

> <mark>FCC Cu is more stable than simple cubic: lower energy</mark>

---

## 🔹 Hexagonal Close-Packed (HCP) Materials

> <mark>FCC & HCP have same packing density, but HCP uses a different lattice</mark>

* Cell vectors:
  $\mathbf{a}_1 = (a, 0, 0), \quad \mathbf{a}_2 = \left( \frac{a}{2}, \frac{a\sqrt{3}}{2}, 0 \right), \quad \mathbf{a}_3 = (0, 0, c)$
* Atoms at:

  * $(0,0,0), (2a/3, 2a/\sqrt{3}, c/2)$

#### Equation 2.7: Fractional Coordinates

$\mathbf{r}_j = \sum_{j=1}^3 f_{ij} \mathbf{a}_j$

> <mark>Useful to describe general atom positions in supercell via lattice vectors.</mark>

### ✅ HCP Geometry

```svg
<svg viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">
  <polygon points="100,50 150,100 100,150 50,100" fill="none" stroke="black"/>
  <circle cx="100" cy="50" r="5" fill="blue"/>
  <circle cx="150" cy="100" r="5" fill="blue"/>
  <circle cx="100" cy="150" r="5" fill="blue"/>
  <circle cx="50" cy="100" r="5" fill="blue"/>
  <circle cx="100" cy="100" r="5" fill="red"/>
</svg>
```

### ✅ HCP Total Energy vs Lattice Ratio (c/a)

**Figure 2.4** – DFT Energy vs $a$ for different $c/a$ ratios:

```svg
<svg viewBox="0 0 300 200" xmlns="http://www.w3.org/2000/svg">
  <line x1="40" y1="10" x2="40" y2="180" stroke="black"/>
  <line x1="40" y1="180" x2="290" y2="180" stroke="black"/>
  <path d="M 50 150 Q 100 80, 150 130 Q 200 100, 250 150" fill="none" stroke="purple" stroke-width="2"/>
  <text x="20" y="20">E</text>
  <text x="250" y="195">c/a</text>
</svg>
```

> Minimum energy at $c/a \approx 1.60$, close to ideal $1.633$. Confirms Cu is not HCP in reality.

---

## 🔹 Crystal Structure Prediction

This section evaluates:

* Whether structures with energy difference $\sim 1 \text{ kJ/mol}$ can be distinguished
* Importance of convergence and controlling all variables (e.g., $c/a$, atom positions)

> <mark>HCP requires varying 2 parameters ($a, c$) vs 1 ($a$) in cubic — more complexity!</mark>

---

## ✅ Summary

* DFT lets us predict equilibrium structure of metals.
* Supercells, primitive cells, and fractional coordinates are essential concepts.
* Graphs show how DFT data fits theoretical models (e.g., Birch–Murnaghan).
* FCC is energetically favored for Cu.

---

## 📌 References

* *Density Functional Theory: A Practical Introduction* – Sholl & Steckel
* Quantum ESPRESSO Docs: [https://www.quantum-espresso.org/](https://www.quantum-espresso.org/)

---
