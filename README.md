# 📊 Quantum ESPRESSO Basics for DFT 

Welcome! This repository is a personal learning journey into using **Quantum ESPRESSO** from scratch for **Density Functional Theory (DFT)** simulations. This README covers everything from basic concepts to your first successful SCF calculation using QE.

---

## 🧠 Overview

### 🔬 What is DFT?

**Density Functional Theory (DFT)** is a quantum mechanical method used to investigate the electronic structure of materials. It simplifies the many-body Schrödinger equation into an equation involving only electron density — making large-scale simulations computationally feasible.

### ⚙️ What is Quantum ESPRESSO?

**Quantum ESPRESSO (QE)** is an open-source suite for electronic structure calculations based on DFT, plane waves, and pseudopotentials. It's widely used for simulating materials at the atomic level.

---

## 📁 Folder Structure

```
.
├── scf.in              # Input file for SCF calculation
├── pseudo/             # Folder containing pseudopotentials
│   └── Si.pz-vbc.UPF
├── tmp/                # Output directory used by QE
├── scf.out             # Output file (generated after running QE)
└── README.md           # This file
```

---

## 💾 SCF Input File Explained (`scf.in`)

Here’s a breakdown of the most important sections in the input file:

### 1️⃣ \&CONTROL

```fortran
&CONTROL
  calculation = 'scf',
  prefix = 'silicon',
  outdir = './tmp',
  pseudo_dir = './pseudo',
  verbosity = 'high',
/
```

Defines calculation type, output directories, and file prefixes.

---

### 2️⃣ \&SYSTEM

```fortran
&SYSTEM
  ibrav = 2,
  celldm(1) = 10.2,
  nat = 2,
  ntyp = 1,
  ecutwfc = 30.0,
  ecutrho = 240.0,
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.01,
/
```

Specifies the Bravais lattice, number of atoms and species, and plane-wave cutoffs.

---

### 3️⃣ \&ELECTRONS

```fortran
&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = 0.7,
/
```

Parameters for the SCF convergence loop.

---

### 4️⃣ ATOMIC\_SPECIES

```bash
ATOMIC_SPECIES
Si  28.0855  Si.pz-vbc.UPF
```

Defines atomic species and associated pseudopotentials.

---

### 5️⃣ ATOMIC\_POSITIONS

```bash
ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
```

Provides atomic coordinates in **crystal** units.

---

### 6️⃣ K\_POINTS

```bash
K_POINTS automatic
6 6 6 0 0 0
```

Defines a Monkhorst-Pack grid for Brillouin zone sampling.

---

## ▶️ How to Run the Simulation

1. Place the `scf.in` file and pseudopotential (`Si.pz-vbc.UPF`) inside the `pseudo/` folder.
2. Open your terminal (WSL on Windows) and run:

```bash
mkdir tmp
pw.x < scf.in > scf.out
```

This will execute a self-consistent DFT calculation for silicon.

---

## 🛠 Requirements

* ✅ Quantum ESPRESSO installed
* ✅ Pseudopotential file in UPF format (`Si.pz-vbc.UPF`)
* ✅ WSL or Linux environment

---

## 📚 References

* [Quantum ESPRESSO Official Site](https://www.quantum-espresso.org/)
* [Quantum ESPRESSO Input Description](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
* [A Primer in Density Functional Theory](https://link.springer.com/book/10.1007/3-540-37072-2)

---

