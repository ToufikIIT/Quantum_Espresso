# 🧮 DFT Calculations for Simple Solids — Full Notes


---

## 🧱 1. Periodic Structures in Solids

**Unit cell** = smallest repeating 3D volume in a solid  
**Lattice vectors** \( \vec{a}_1, \vec{a}_2, \vec{a}_3 \) describe periodicity:

- Simple cubic: \( \vec{r}_{n_1, n_2, n_3} = n_1 a \hat{x} + n_2 a \hat{y} + n_3 a \hat{z} \)
- **Supercell**: Large collection of unit cells used in DFT

Atoms repeat using periodic boundary conditions (PBC), making DFT simulations feasible.

---

## 📉 2. Lattice Constant and Energy Optimization

To find equilibrium structure:
1. Vary lattice parameter \( a \)
2. Compute total energy \( E_{tot}(a) \)
3. Fit to a smooth curve

### ➕ Taylor Expansion:
\[
E_{tot}(a) \approx E_0 + \alpha(a - a_0) + \beta(a - a_0)^2
\]

Simplified fit:
\[
E_{tot}(a) = E_0 + \beta(a - a_0)^2
\]

Where:
- \( a_0 \): Optimal lattice constant
- \( B_0 \): Bulk modulus
- \( E_0 \): Minimum energy

### ⚙️ Birch–Murnaghan Equation:
Used for more accurate fitting:

\[
E_{tot}(a) = E_0 + \frac{9 V_0 B_0}{16} \left\{ \left[\left(\frac{a_0}{a}\right)^2 - 1\right]^3 B'_0 + \left[\left(\frac{a_0}{a}\right)^2 - 1\right]^2 \left[6 - 4\left(\frac{a_0}{a}\right)^2\right] \right\}
\]

---

## 💠 3. FCC Structure (Face-Centered Cubic)

Common for Cu and other metals.

- Atom positions in conventional cell:
  \[
  (0, 0, 0),\quad (a/2, a/2, 0),\quad (a/2, 0, a/2),\quad (0, a/2, a/2)
  \]

- Primitive lattice vectors:
  \[
  \vec{a}_1 = a \left(\frac{1}{2}, 0, \frac{1}{2}\right),\quad
  \vec{a}_2 = a \left(0, \frac{1}{2}, \frac{1}{2}\right),\quad
  \vec{a}_3 = a \left(\frac{1}{2}, \frac{1}{2}, 0\right)
  \]

- Nearest neighbor distance: \( a/\sqrt{2} \)

### ✅ FCC Cu DFT Results:
- Energy minimized at \( a = 3.64 \) Å
- Bulk modulus: \( B_0 = 142 \) GPa

---

## ⬢ 4. HCP Structure (Hexagonal Close-Packed)

HCP and FCC have same packing efficiency (74%) but different geometry.

- Lattice vectors:
  \[
  \vec{a}_1 = (a, 0, 0),\quad
  \vec{a}_2 = \left(\frac{a}{2}, \frac{a \sqrt{3}}{2}, 0\right),\quad
  \vec{a}_3 = (0, 0, c)
  \]

- Atom positions:
  \[
  (0,0,0),\quad \left(\frac{2a}{3}, \frac{2a}{\sqrt{3}}, \frac{c}{2}\right)
  \]

- **Fractional coordinates**: 
  \[
  \vec{r}_i = \sum_{j=1}^3 f_{ij} \vec{a}_j
  \]

### 🔍 c/a Ratio Effect:
- Ideal \( c/a = 1.633 \)
- DFT finds lowest energy near \( c/a = 1.60 \)
- Plot of \( E_{tot}(a) \) for several \( c/a \) ratios used

---

## 🧠 5. Crystal Structure Prediction Strategy

To determine which structure (FCC or HCP) is more stable:

- Use DFT to compute total energy for each structure
- Compare energies at optimized lattice constants
- Cu prefers **FCC** due to slightly lower energy

> ⚠️ Energy differences are small (1–5 kJ/mol), so precision is critical.


---

## 📚 References

- Sholl & Steckel, *Density Functional Theory: A Practical Introduction*
- Quantum ESPRESSO Project: https://www.quantum-espresso.org/

---

