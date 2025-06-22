
# 📘Reciprocal Space, *k*-Points & Pseudopotentials

This chapter explores essential numerical foundations for periodic DFT: reciprocal space, Brillouin zones, k-point sampling, energy cutoffs, and pseudopotentials. Each section includes intuitive explanations, annotated equations, and best practices for real-world simulations.

---

## 🧭 1.1 Reciprocal Space and *k*-Points

### Why Reciprocal Space?

Solids are periodic. Real space isn’t the most natural setting to describe periodicity — **reciprocal space** is. Many physical quantities, such as electronic band structure, are most naturally expressed in terms of reciprocal vectors and wavevectors $\\vec{k}$
.

> 🔁 Reciprocal space simplifies operations like solving the Schrödinger equation, integrating over the unit cell, and describing interference patterns.

---

## 📐 1.1.1 Plane Waves and the Brillouin Zone

### Bloch’s Theorem

In a periodic potential:

$$
\phi_k(\vec{r}) = e^{i\vec{k} \cdot \vec{r}} u_k(\vec{r})
$$

- $\\ u_k(\vec{r}) $: periodic with the lattice
- $\\ \vec{k} $: crystal momentum

This allows the Schrödinger equation to be solved for each \( \vec{k} \) separately.

---

### Reciprocal Lattice Vectors

For primitive vectors \( \vec{a}_1, \vec{a}_2, \vec{a}_3 \), define reciprocal vectors:

$$
\vec{b}_1 = 2\pi \frac{\vec{a}_2 \times \vec{a}_3}{\vec{a}_1 \cdot (\vec{a}_2 \times \vec{a}_3)} \\
\vec{b}_2 = 2\pi \frac{\vec{a}_3 \times \vec{a}_1}{\vec{a}_1 \cdot (\vec{a}_2 \times \vec{a}_3)} \\
\vec{b}_3 = 2\pi \frac{\vec{a}_1 \times \vec{a}_2}{\vec{a}_1 \cdot (\vec{a}_2 \times \vec{a}_3)}
$$

These satisfy:

![Equation](https://latex.codecogs.com/png?\fg{ffffff}\vec{a}_i%20\cdot%20\vec{b}_j%20=%202\pi%20\delta_{ij})


> 🧠 These relations define the dual lattice used in wavevector-based descriptions.

---

### Brillouin Zone

The **first Brillouin zone** (BZ) is the Wigner-Seitz cell in reciprocal space.

**Volume:**

$$
V_{BZ} = \frac{(2\pi)^3}{V_{cell}}
$$

---

## 📊 1.1.2 Integration in *k*-Space

DFT calculations need integrals over the Brillouin zone:

$$
g = \frac{V_{cell}}{(2\pi)^3} \int_{BZ} g(\vec{k}) \, d\vec{k}
$$

These are approximated by **sampling *k*-points** on a grid.

---

### Integration Techniques

#### Trapezoidal Rule:

$$
\int_{-1}^{1} f(x) dx \approx \frac{h}{2} [f(-1) + f(1) + 2\sum f(x_j)]
$$

#### Gaussian Quadrature:

$$
\int_{-1}^{1} f(x) dx \approx \sum c_j f(x_j)
$$

> 🔍 Useful for accurate integration with fewer points

---

## 🧮 1.1.3 Monkhorst–Pack Grids

Instead of random sampling, we use systematic *k*-point meshes.

### Monkhorst–Pack Grids:

- Construct uniform grids $\\ M \times M \times M $
- Exploit symmetry → sample **irreducible Brillouin zone (IBZ)**

#### Table: Total Energy vs Grid Density

| M  | Energy (eV) | IBZ Points | Time Ratio |
|----|-------------|------------|------------|
| 1  | -3.8061     | 1          | 1.0        |
| 8  | -3.7661     | 56         | 31.2       |
| 14 | -3.7659     | 84         | 39.7       |

> ✅ Convergence achieved by M = 8  
> 🔁 Odd M includes Γ-point

---

## 🔪 1.2 Energy Cutoff for Plane Waves

Plane wave basis:

Our lengthy discussion of k space began with Bloch’s theorem, which tells us that solutions of the Schro¨dinger equation for a supercell have the form

$$
\phi_k(\vec{r}) = e^{i\vec{k} \cdot \vec{r}} u_k(\vec{r})
$$

where uk(r) is periodic in space with the same periodicity as the supercell.

It is now time to look at this part of the problem more carefully. The periodicity of uk(r) means that it can be expanded in terms of a special set of plane waves

![Equation](https://latex.codecogs.com/png?\fg{ffffff}u_k(\mathbf{r})%20=%20\sum_{\mathbf{G}}%20c_{\mathbf{G}}%20\exp[i%20\mathbf{G}%20\cdot%20\mathbf{r}])

 where the summationis overall vectors defined by ![Equation](https://latex.codecogs.com/png?\fg{ffffff}\mathbf{G}%20=%20m_1\mathbf{b}_1%20+%20m_2\mathbf{b}_2%20+%20m_3\mathbf{b}_3)
 with integer values for mi. These set of vectors defined by G in reciprocal
 space are defined so that for any real space lattice vector ai, ![Equation](https://latex.codecogs.com/png?\fg{ffffff}\mathbf{G}\cdot\mathbf{a}_i%20=%202\pi%20m)

 Combining the two equations above gives

$$
\phi_k(\vec{r}) = \sum_{\vec{G}} c_{\vec{k} + \vec{G}} e^{i(\vec{k} + \vec{G}) \cdot \vec{r}}
$$

To limit the infinite sum, use:

### Cutoff Condition:

$$
\frac{\hbar^2}{2m} |\vec{k} + \vec{G}|^2 < E_{cut}
$$

Or equivalently:

$$
|\vec{k} + \vec{G}|^2 < G_{cut}^2
$$

> ⚠️ Low cutoffs lead to inaccuracy  
> ✅ Higher cutoffs = more precision, more cost

---

## ⚛️ 1.3 Pseudopotentials

### Why Pseudopotentials?

Core electron wavefunctions are highly oscillatory → need **many plane waves**. But:

> ⚠️ Core electrons don't strongly influence chemical bonding.

Hence, we **approximate core effects** with pseudopotentials → smooth wavefunctions → **lower cutoff**.

---

### Frozen Core Approximation

- Replace core electrons with fixed smoothed density
- Only valence electrons are explicitly treated
- Used in nearly all plane-wave DFT codes

**Opposite**: All-electron calculations (rare in PW codes)

---

### Transferability

A good pseudopotential:
- Built from an isolated atom
- Re-usable in any compound
- Needs **no tuning per compound**

> 🔁 This ability is called **transferability**

---

### Hard vs. Soft Pseudopotentials

| Type | Description                        | Cutoff Needed |
|------|------------------------------------|----------------|
| Hard | Highly accurate but oscillatory    | High           |
| Soft | Smoother, easier to converge       | Low            |

---

### Ultrasoft Pseudopotentials (USPP)

- Developed by Vanderbilt
- Use **fewer plane waves**
- Require **empirical tuning**
- Available for most elements in DFT libraries

> ✅ Efficient but must be tested carefully

---

### Projector Augmented-Wave (PAW)

- Invented by Blöchl
- Combines **all-electron accuracy** + **plane-wave efficiency**
- More accurate than USPP
- Used in VASP and QE

---

### Method Comparison Table

| Method         | Core Electrons | Accuracy      | Efficiency | Notes                                |
|----------------|----------------|---------------|------------|--------------------------------------|
| All-Electron   | Explicit       | 🔵 Very High  | 🔴 Slow     | Full description of atom             |
| USPP           | Frozen         | 🟡 Medium     | 🟢 Fast     | May require tuning                   |
| PAW            | Frozen (augmented) | 🔵 High | 🟢 Fast     | Best mix of speed & accuracy         |

---

## 🧠 Final Summary

- Reciprocal space describes periodic systems better
- Brillouin Zone integrals are done via discrete *k*-points
- Plane wave expansions are truncated using an energy cutoff
- Pseudopotentials make calculations feasible for large systems

---

## ✅ Best Practices

- 🔍 Always test **convergence** wrt *k*-points and cutoff
- 🧪 Use **odd Monkhorst–Pack** grids to include Γ-point
- ⚖️ Balance **soft pseudopotentials** with accuracy needs
- 🚀 Prefer **PAW** for accuracy + speed

---
