# MBD Framework: Many-Body Dispersion Computations

**`mbd-framework`** is a computational toolkit that calculates precise Many-Body Dispersion (MBD) scaling bounds and Tkatchenko-Scheffler (TS) screening parameters. It acts as the numerical bridge between formal mathematical proofs (derived in Lean 4) and first-principles quantum chemistry computations (powered by PySCF).

---

## 🚀 Installation

Install the framework globally via PyPI:

```bash
pip install mbd-framework
```

---

## 💻 Quick Start & Command Line Tools

Once installed, the framework provides three native command-line tools to instantly run calculations on your molecular targets:

### 1. Extract Atomic Density Bounds (`mbd-compute`)
Run background PySCF calculations to extract absolute finite-field Cartesian dipole tensors, compute precise atomic polarizabilities ($\alpha$), and map them to the universal TS scaling parameter $x = V_{\text{Bohr}} / \alpha$. The results are safely checkpointed into a local `database.json`.

```bash
# Example: Extract the bounds for Benzene using the aug-cc-pVDZ basis set
mbd-compute --molecule Benzene --basis aug-cc-pVDZ

*Note: Cross-validation is on by default. Use `--methods rhf` for HF only, or `--methods rks-pbe,rks-b3lyp` to compare functionals.*

*Supported molecules:* `Benzene`, `Naphthalene`, `He`, `Ne`, `Xe`

### 2. Simulate Crystal Dispersion (`mbd-crystal`)
Generate a structurally symmetric `7x7x7` Cartesian atomic lattice array (thousands of interacting pairs) to compute the macroscopic isotropic dispersion grid ($C_{6} \cdot \varepsilon^{-x} / R^{6}$) against your extracted bounds, yielding total empirical lattice energies in `kJ/mol`.

```bash
# Example: Simulate Benzene dispersion scaling in a vacuum (epsilon = 1.0)
mbd-crystal --target Benzene --epsilon 1.0
```

### 3. SERS Mathematical Equivalence (`mbd-sers`)
A strict numerical checking endpoint that compares the macroscopic structural Surface-Enhanced Raman Scattering (SERS) exponential quenching envelope ($\exp(-\rho)$) against the framework's intrinsic MBD interacting boundaries ($\varepsilon^{-x}$).

```bash
# Example: Test SERS equivalence for Benzene
mbd-sers --target Benzene --epsilon 2.0
```

---

## 🛠 Features at a Glance
- **Strictly Typed Physics**: Detects and corrects SCF convergence to saddle points via Pople stability iteration before computing properties.
- **Automated Checkpointing**: Computations automatically stack results into a local `database.json` to prevent losing time-intensive finite-field SCF runs.
- **Seamless Academic Pipeline**: Numerically validates mathematical limit theorems bounding dispersion properties for large-scale Molecular Crystals.

---

*For detailed physics explanations, API structures, and Lean 4 formal mathematical proofs, please visit the main academic repository at [GitHub: MBD-Theoretical-Framework](https://github.com/edqa/MBD-Theoretical-Framework).*
