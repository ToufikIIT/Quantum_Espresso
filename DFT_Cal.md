# 📘 DFT Calculations for Simple Solids

---

## 📌 Table of Contents

1. [Introduction](#introduction)
2. [Understanding Unit Cells and Lattice Parameters](#understanding-unit-cells-and-lattice-parameters)
3. [Total Energy and Lattice Constant: How to Find the Best Structure](#total-energy-and-lattice-constant-how-to-find-the-best-structure)
4. [FCC Structure: A Case Study with Copper](#fcc-structure-a-case-study-with-copper)
5. [HCP Structure: Another Common Packing Type](#hcp-structure-another-common-packing-type)


---

## 🔷 Introduction

**Why DFT?**
Density Functional Theory (DFT) allows us to compute the **total energy** of a material based on how atoms are arranged in it. This is useful to predict the **most stable structure** of solids — the one with the **lowest energy**.

Think of it as testing different possible structures (like stacking blocks in different ways) and calculating which arrangement is most natural (energetically preferred).

---

## 🔹 Understanding Unit Cells and Lattice Parameters

### What is a Unit Cell?

It is the smallest 3D repeating block of a crystal. When repeated in space, it forms the entire structure of the solid.

### Simple Cubic Example:

* **Atoms at corners** only.
* Lattice parameter (edge length) = $a$
* Any atomic position can be described by:

$$
\vec{r} = (n_1 a, n_2 a, n_3 a) \text{ where } n_i \in \mathbb{Z}
$$

A **supercell** is a large group of these unit cells, used in DFT to simulate bulk materials more accurately.

---

## 🔹 Total Energy and Lattice Constant: How to Find the Best Structure

We want to find the value of $a$ (lattice constant) that makes the total energy $E_{tot}$ of the material **minimum**.

### Step-by-step:

1. Choose different values of $a$
2. For each $a$, calculate $E_{tot}$ using DFT
3. Fit a curve to the points
4. Find the $a$ value where energy is lowest → that’s the best lattice constant

### 🔸 Taylor Expansion (Simple Fit)

e.g., near the minimum:

$$
E_{tot}(a) = E_0 + \beta(a - a_0)^2
$$

Where:

* $a_0$: equilibrium lattice constant (e.g., 2.43 Å)
* $\beta$: curvature; related to how stiff the lattice is

### 🔸 Birch-Murnaghan Equation (Advanced Fit)

This is used for better accuracy near the minimum:

$$
E_{\text{tot}}(a) = E_0 + \frac{9V_0 B_0}{16} \left( \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^3 B_0' + \left[ \left( \frac{a_0}{a} \right)^2 - 1 \right]^2 \left[ 6 - 4 \left( \frac{a_0}{a} \right)^2 \right] \right)
$$

Where:

* $E_0$: energy at equilibrium
* $V_0$: volume at equilibrium
* $B_0$: bulk modulus (resistance to compression)
* $B_0'$: pressure derivative of bulk modulus

---

## 🔹 FCC Structure: A Case Study with Copper

### FCC (Face-Centered Cubic):

* Atoms located at corners and face centers of the cube

**Coordinates of atoms:**

* (0, 0, 0)
* ($a/2, a/2, 0$)
* ($a/2, 0, a/2$)
* ($0, a/2, a/2$)

**Primitive Vectors:**

$$
\vec{a}_1 = a\left(\frac{1}{2}, 0, \frac{1}{2}\right),\quad
\vec{a}_2 = a\left(0, \frac{1}{2}, \frac{1}{2}\right),\quad
\vec{a}_3 = a\left(\frac{1}{2}, \frac{1}{2}, 0\right)
$$

### Energy Calculation:

* The curve of $E_{tot}$ vs. $a$ shows a minimum at $a = 3.64 \text{ Å}$
* This means 3.64 Å is the most stable size for the FCC cell of Cu
* Bulk modulus from this fit: $142 \text{ GPa}$

✅ **Conclusion:** FCC is the stable structure for copper.

---

## 🔹 HCP Structure: Another Common Packing Type

**Hexagonal Close-Packed (HCP)** is another way to stack atoms compactly:

* 2 lattice parameters: $a$ and $c$
* Ideal ratio $c/a \approx 1.633$

**Primitive Vectors:**

$$
\vec{a}_1 = (a, 0, 0), \quad
\vec{a}_2 = \left(\frac{a}{2}, \frac{a\sqrt{3}}{2}, 0\right), \quad
\vec{a}_3 = (0, 0, c)
$$

**Atom Positions:**

* (0, 0, 0)
* $\left(\frac{2a}{3}, \frac{2a}{\sqrt{3}}, \frac{c}{2}\right)$

### Fractional Coordinates

Fractional positions help describe atomic locations within arbitrary unit cells:

$$
\[
\vec{r}_j = \sum_{i=1}^3 f_{ij} \vec{a}_i
\]
$$


Where $f_{ij} \in [0, 1)$

### HCP Energy vs. $c/a$

* For each $c/a$ ratio, calculate total energy
* Curve shows minimum near $c/a = 1.60$, close to ideal value
* Indicates that HCP is a **competitive** structure for some materials

---

## 🔹 Comparing Structures to Predict Stability

We calculate total energy for each structure and compare:

| Structure | Parameters            | Energy Comparison             |
| --------- | --------------------- | ----------------------------- |
| FCC       | 1 parameter ($a$)     | Lower energy → stable         |
| HCP       | 2 parameters ($a, c$) | Slightly higher energy for Cu |

> ⚠️ Even small energy differences (few kJ/mol) matter — DFT must be **highly accurate**

Convergence tests are important:

* Fine k-point mesh
* Proper energy cutoff
* Sufficient number of atoms

---

## 📖 References

* Sholl & Steckel, *Density Functional Theory: A Practical Introduction*
* Quantum ESPRESSO: [https://www.quantum-espresso.org](https://www.quantum-espresso.org)

---
