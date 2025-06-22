
# 🧪 Quantum ESPRESSO Input File: Detailed Explanation

Quantum ESPRESSO (QE) uses plain-text input files to define **simulation parameters, system configuration, and control options** for running calculations. Here's a complete breakdown of an input file used with `pw.x` (Plane-Wave Self-Consistent Field code).

---

## 🔹 General Structure

A typical QE input file is divided into the following parts:

```fortran
&CONTROL
/
&SYSTEM
/
&ELECTRONS
/
&IONS          (optional)
/
&CELL          (optional)
/ 

ATOMIC_SPECIES
ATOMIC_POSITIONS {crystal|alat|bohr|angstrom}
K_POINTS {automatic|gamma|crystal}
CELL_PARAMETERS {alat|bohr|angstrom}   (optional)
```

---

## ✅ &CONTROL Section

Controls the overall run — defines what kind of calculation is performed, where files are stored, and how much detail is printed. Below are the key variables and their purposes:

- `calculation`: Type of task to be performed.
  - Options: `'scf'`, `'nscf'`, `'bands'`, `'relax'`, `'md'`, `'vc-relax'`, `'vc-md'`
  - **Default**: `'scf'`

- `prefix`: Name prefix for output files.

- `outdir`: Temporary directory for output.
  - **Default**: Uses the `ESPRESSO_TMPDIR` environment variable if set, otherwise `./`.

- `pseudo_dir`: Where pseudopotentials are stored.

- `verbosity`: Level of detail in output.
  - Options: `'low'`, `'high'`, `'default'`, `'minimal'`, `'debug'`
  - `'debug'` and `'medium'` → `'high'`
  - `'default'` and `'minimal'` → `'low'`

- `wf_collect`: Logical. **Obsolete** — no longer implemented.

- `title`: Title printed in output.
  - **Default**: `''`

- `restart_mode`: How to restart a calculation.
  - Options: `'from_scratch'` (default), `'restart'`
  - `'restart'` resumes an interrupted calculation.

- `nstep`: Number of MD or relaxation steps.
  - **Default**: `1` for `'scf'`, `'nscf'`, `'bands'`; `50` otherwise
  - Set `nstep = 0` for a "dry run" (initialization only).

- `disk_io`: Disk I/O control.
  - `'high'`: Save everything
  - `'medium'`: Save charge and selectively save WFs
  - `'low'`: Save WFs only at the end (default for `'scf'`)
  - `'nowf'`, `'minimal'`, `'none'`: Save progressively less

| Keyword       | Meaning                                                                 |
|---------------|-------------------------------------------------------------------------|
| `calculation` | Type of calculation: `'scf'`, `'relax'`, `'vc-relax'`, `'nscf'`, etc.   |
| `prefix`      | Name prefix for output files (e.g., `prefix.save`).                     |
| `pseudo_dir`  | Directory containing pseudopotential files.                             |
| `outdir`      | Temporary directory for wavefunctions, charge density, etc.             |
| `tstress`     | Calculate stress tensor? (`.true.` or `.false.`)                        |
| `tprnfor`     | Print forces? (`.true.` or `.false.`)                                   |
| `verbosity`   | Level of output detail: `'low'`, `'default'`, `'high'`, `'debug'`       |
| `wf_collect`  | Collect wavefunctions into a single file (important for parallel runs)  |

**Example:**
```fortran
&CONTROL
  calculation = 'scf',
  prefix = 'silicon',
  outdir = './tmp/',
  pseudo_dir = './pseudo/',
  tstress = .true.,
  tprnfor = .true.,
  wf_collect = .true.
/
```

---

## ⚙️ &SYSTEM Section

Defines the physical properties of the system, including lattice, number of atoms, cutoff energies, and smearing settings:

- `ibrav`: Bravais lattice index.
  - `ibrav = 0`: Use `CELL_PARAMETERS` instead of automatic lattice
  - Use either `celldm(i)` or `{A, B, C, cosAB, cosAC, cosBC}`

- `celldm(i)`: Lattice parameters in Bohr units.
  - `celldm(1)`: Lattice constant
  - Others: Ratios and angles depending on `ibrav`

- `A, B, C, cosAB, cosAC, cosBC`: Optional lattice constants in Angstrom (used when `ibrav ≠ 0`).

- `nat`: Total number of atoms in the unit cell.

- `ntyp`: Number of distinct atomic species.

- `tot_charge`: Total charge of the cell.
  - `0`: Neutral (default)
  - `+1`: Missing 1 electron
  - `-1`: Extra 1 electron

- `ecutwfc`: Plane-wave cutoff for wavefunctions (in Ry).

- `ecutrho`: Cutoff for charge density. Typically 4–8× `ecutwfc`.

- `occupations`: How electron occupancy is assigned.
  - `'fixed'`: For insulators
  - `'smearing'`: For metals
  - `'tetrahedra'`: For improved BZ integration

- `smearing`: Smearing method if applicable.
  - Examples: `'gaussian'`, `'methfessel-paxton'`, `'fermi-dirac'`

- `degauss`: Smearing width in Ry.

| Keyword                  | Description                                                      |
|--------------------------|------------------------------------------------------------------|
| `ibrav`                  | Bravais lattice index (0 for custom lattice)                     |
| `celldm(i)`              | Lattice parameters (used with `ibrav > 0`)                        |
| `nat`                    | Number of atoms in the unit cell                                 |
| `ntyp`                   | Number of atomic species                                         |
| `ecutwfc`                | Plane-wave cutoff energy (in Ry)                                 |
| `ecutrho`                | Cutoff for charge density and potential (usually ≥ 4×ecutwfc)    |
| `occupations`            | Type of occupation: `'smearing'`, `'fixed'`, `'tetrahedra'`      |
| `smearing`               | Smearing method: `'gaussian'`, `'mp'`, `'fermi-dirac'`, etc.     |
| `degauss`                | Smearing width (in Ry)                                           |
| `nspin`                  | 1 for non-magnetic, 2 for spin-polarized                         |
| `starting_magnetization`| Initial spin polarization guess                                   |

**Example:**
```fortran
&SYSTEM
  ibrav = 2,
  celldm(1) = 10.2,
  nat = 2,
  ntyp = 1,
  ecutwfc = 40.0,
  ecutrho = 320.0,
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.01
/
```

---

## ⚡ &ELECTRONS Section

Controls **SCF (Self-Consistent Field) convergence and mixing**.

| Keyword           | Description                                                    |
|-------------------|----------------------------------------------------------------|
| `conv_thr`        | Convergence threshold for SCF loop (in Ry)                     |
| `mixing_beta`     | Mixing factor for charge density (0.1–0.7 typical)             |
| `electron_maxstep`| Maximum number of SCF steps                                    |
| `diagonalization` | Diagonalization method (e.g., `'david'`, `'cg'`, `'rmm-davidson'`) |

**Example:**
```fortran
&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = 0.7
/
```

---

## 🧲 &IONS Section *(optional)*

Controls **ionic relaxation**, used in `'relax'` or `'vc-relax'` calculations.

| Keyword         | Description                         |
|------------------|-------------------------------------|
| `ion_dynamics`   | Method: `'bfgs'`, `'damp'`, `'md'`  |

```fortran
&IONS
  ion_dynamics = 'bfgs'
/
```

---

## 🧱 &CELL Section *(optional)*

Defines **variable-cell optimization** settings.

| Keyword         | Description                         |
|------------------|-------------------------------------|
| `cell_dynamics`  | `'bfgs'`, `'damp'`, `'pr'`, `'md'` |
| `press`          | External pressure (in kbar)         |
| `cell_factor`    | Scaling factor for cell variation   |

```fortran
&CELL
  cell_dynamics = 'bfgs',
  press = 0.0,
  cell_factor = 2.0
/
```

---

## 🧪 ATOMIC_SPECIES

Defines each **element type** used, its **atomic mass**, and **pseudopotential file**.

```fortran
ATOMIC_SPECIES
Si 28.0855 Si.pz-vbc.UPF
```

---

## ⚛️ ATOMIC_POSITIONS

Specifies **atomic coordinates** and format.

- Coordinate types: `{crystal}`, `{alat}`, `{bohr}`, `{angstrom}`

```fortran
ATOMIC_POSITIONS crystal
Si  0.00  0.00  0.00
Si  0.25  0.25  0.25
```

---

## 🧾 K_POINTS

Defines the **k-point grid** for Brillouin zone integration.

### Option 1: Automatic Grid
```fortran
K_POINTS automatic
4 4 4 0 0 0
```

### Option 2: Gamma-only Point
```fortran
K_POINTS gamma
```

---

## 📐 CELL_PARAMETERS *(if `ibrav = 0`)*

Provides custom **lattice vectors**.

```fortran
CELL_PARAMETERS angstrom
5.43 0.00 0.00
0.00 5.43 0.00
0.00 0.00 5.43
```

---

## 🧩 Complete Example: SCF Input for Silicon

```fortran
&CONTROL
  calculation = 'scf',
  prefix = 'silicon',
  outdir = './tmp/',
  pseudo_dir = './pseudo/',
  tstress = .true.,
  tprnfor = .true.,
  wf_collect = .true.
/

&SYSTEM
  ibrav = 2,
  celldm(1) = 10.2,
  nat = 2,
  ntyp = 1,
  ecutwfc = 40.0,
  ecutrho = 320.0,
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.01
/

&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = 0.7
/

ATOMIC_SPECIES
Si 28.0855 Si.pz-vbc.UPF

ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25

K_POINTS automatic
4 4 4 0 0 0
```

---

## 📚 References

- [QE Official Documentation](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
- [QE Input Generator](https://www.materialscloud.org/work/tools/qeinputgenerator)
- [Pseudopotentials Library](https://www.quantum-espresso.org/pseudopotentials)

---

## 📘 Additional Details from Official Documentation

### 🔹 &CONTROL Section (Expanded)

Controls the overall run: what calculation is performed, where files are saved, etc.

- `calculation`: Type of task.
  - Options: `'scf'`, `'nscf'`, `'bands'`, `'relax'`, `'md'`, `'vc-relax'`, `'vc-md'` (Default: `'scf'`)

- `prefix`: Prefix for output files.

- `outdir`: Temporary directory for files.
  - Default: `ESPRESSO_TMPDIR` env variable if set, else `./`.

- `pseudo_dir`: Directory of pseudopotentials.

- `verbosity`: Level of output.
  - Options: `'high'`, `'low'`, `'default'`, `'minimal'`, `'debug'`

- `wf_collect`: *Obsolete* — no longer used.

- `title`: String printed in output.

- `restart_mode`: `'from_scratch'` (default), `'restart'`

- `nstep`: Number of MD or optimization steps (Default: 1 for `scf`, `nscf`; 50 otherwise)

- `disk_io`: Controls disk usage during the run.
  - `'high'`: Max disk usage; save everything.
  - `'medium'`: Save charge and WFs depending on k-points.
  - `'low'`: Save WFs only at the end; default for `scf`.
  - `'nowf'`, `'minimal'`, `'none'`: Progressive levels of saving less to disk.

---

### 🔹 &SYSTEM Section (Expanded)

- `ibrav`: Bravais-lattice index.
  - If `ibrav = 0`, use `CELL_PARAMETERS`.
  - Use either `celldm(i)` or `{A, B, C, cosAB, cosAC, cosBC}`.

- `celldm(i)`: Lattice parameters in Bohr.

- `A, B, C, cosAB, cosAC, cosBC`: Used for defining lattice in Angstrom when `ibrav ≠ 0`.

- `nat`: Number of atoms in the unit cell.

- `ntyp`: Number of atomic species.

- `tot_charge`: Total charge (e.g., `+1` = remove 1e⁻, `-1` = add 1e⁻). Default: 0 (neutral).

- `ecutwfc`: Cutoff energy for wavefunctions (in Ry).

- `ecutrho`: Cutoff energy for charge density. Typically 4–8× `ecutwfc`.

- `occupations`: Specifies how electrons are filled.
  - Options: `'fixed'`, `'smearing'`, `'tetrahedra'`

- `smearing`: Method of smearing (e.g., `'gaussian'`, `'mp'`, `'fermi-dirac'`).

- `degauss`: Width of smearing (in Ry). Balance convergence and accuracy.

---
