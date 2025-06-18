# 📘 DFT Calculations for Simple Solids 

---

## 📌 Contents

1. [Introduction](#introduction)
2. [Periodic Structures and Lattice Parameters](#periodic-structures-and-lattice-parameters)
3. [Total Energy and Lattice Constant](#total-energy-and-lattice-constant)
4. [FCC Structure](#fcc-structure)
5. [HCP Structure](#hcp-structure)
---

## 🔷 Introduction

DFT allows us to predict the **equilibrium crystal structure** of pure solid materials (e.g., metals). It focuses on computing the **total energy** of a material as a function of its atomic arrangement.

We want to determine the **positions of all atoms** within a repeating 3D crystal using concepts like **unit cells**, **supercells**, and **lattice vectors**.

---

## 🔹 Periodic Structures and Lattice Parameters

A **simple cubic structure**:

* Unit cell: cube with one atom at each corner
* Edge length: $a$
* Atomic positions: $\vec{r} = (n_1a, n_2a, n_3a)$

The **supercell** is a 3D volume made of these repeated unit cells, used for DFT calculations.

---

## 🔹 Total Energy and Lattice Constant

To find the most stable crystal structure, we vary the lattice constant $a$ and compute the total energy $E_{tot}$:

### 📐 Taylor Expansion:

$$
E_{tot}(a) \approx E_{tot}(a_0) + \alpha (a - a_0) + \beta (a - a_0)^2
$$

Often simplified by fitting:

$$
E_{tot}(a) = E_0 + \beta(a - a_0)^2
$$

Where:

* $a_0 = 2.43 \text{ Å}$
* $B_0 = 103 \text{ GPa}$ (bulk modulus)

### 📐 Birch-Murnaghan Equation:

$$
E_{\text{tot}}(a) = E_0 + \frac{9V_0 B_0}{16} \left( \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^3 B_0' + \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^2 \left[ 6 - 4 \left( \frac{a_0}{a} \right)^2 \right] \right)
$$


This allows broader and more accurate fitting of DFT data.

---

## 🔹 FCC Structure

**Face-Centered Cubic (FCC)** is common for metals like Cu:

* Atom positions: $(0,0,0), (a/2,a/2,0), (a/2,0,a/2), (0,a/2,a/2)$
* Primitive vectors:

$$
\vec{a}_1 = a \left(\frac{1}{2}, 0, \frac{1}{2}\right),\quad \vec{a}_2 = a \left(0, \frac{1}{2}, \frac{1}{2}\right),\quad \vec{a}_3 = a \left(\frac{1}{2}, \frac{1}{2}, 0\right)
$$

* Nearest neighbor distance: $a/\sqrt{2}$

### ⚙️ FCC Lattice Energy

* Minimum total energy at $a = 3.64 \text{ Å}$
* Bulk modulus $B_0 = 142 \text{ GPa}$

> FCC Cu is energetically more stable than simple cubic Cu.

---

## 🔹 HCP Structure

**Hexagonal Close-Packed (HCP)** shares same packing efficiency as FCC, but differs geometrically:

* Cell vectors:

$$
\vec{a}_1 = (a, 0, 0),\quad \vec{a}_2 = \left(\frac{a}{2}, \frac{a\sqrt{3}}{2}, 0\right),\quad \vec{a}_3 = (0, 0, c)
$$

* Atom positions:

$$
(0,0,0), \quad \left(\frac{2a}{3}, \frac{2a}{\sqrt{3}}, \frac{c}{2}\right)
$$

### 🔹 Fractional Coordinates:

$$
\vec{r}_j = \sum_{j=1}^3 f_{ij} \vec{a}_j
$$

Useful for placing atoms inside complex cells using scaled lattice vectors.

### ⚙️ HCP Energy vs c/a

* DFT total energy plotted vs. $a$, for various $c/a$ values
* Energy minimized around $c/a = 1.60$, close to ideal value 1.633

---

## 📖 References

* Sholl & Steckel, *DFT: A Practical Introduction*
* Quantum ESPRESSO official docs: [https://www.quantum-espresso.org/](https://www.quantum-espresso.org/)

---
