# Formal Verification of a Specific Many-Body Dispersion (MBD) Bound

## Problem
While Many-Body Dispersion (MBD) and Tkatchenko-Scheffler (TS) screening methods are heavily utilized in computational chemistry, their core mathematical screening bounds often lack formal, machine-checked verification. As dispersion corrections increasingly dictate the accuracy of massive molecular crystal and material simulations, ensuring the absolute mathematical rigor of the underlying scaling limits (e.g., $x = V_{\text{Bohr}} / \alpha$) is critical.

## Existing Methods
Large-scale numerical integration of Many-Body Dispersion is currently handled efficiently by highly optimized, established libraries such as [libMBD](https://github.com/libmbd/libmbd) and PyMBD. These tools are the gold standard for producing dispersion energies across thousands of interacting atomic pairs in realistic materials. 

## Our Contribution
We are not replacing operational implementations like `libMBD`. Instead, our contribution is **strictly scoped** to providing the first *machine-checked Lean 4 formalization* of the specific $x = V_{\text{Bohr}} / \alpha$ proportionality limit. 

This repository unites **interactive theorem proving** with **first-principles quantum chemistry**:
1. **Lean 4 Proofs**: We mathematically prove that MBD screening establishes a formal parallel to the SERS exponential quenching envelope $\exp(-\rho)$, and verify the bounding theorem $\varepsilon^{-x} \le 1$.
2. **Strict Python Extraction**: We provide a Python computational bridge (`mbd-framework`) that interfaces with PySCF to extract absolute finite-field Cartesian dipole tensors, strictly enforcing the Lean-verified bounds on the resulting physical parameters.

## Results and Validation
Our strictly typed Python bridge successfully bounds atomic polarizabilities against the universal Bohr volume ($V_{\text{Bohr}} \approx 0.62 \text{ \AA}^3$).

### A. Baseline $x$ Ratios (PySCF Validation)
Using `aug-cc-pVDZ` configurations, we validated the extraction against expected TS empirical standards. Strongly-bound electron systems assert heavier screening constraints, while diffuse atoms assert minimal screening:
- **Helium ($He$)**: $x = 3.26$ 
- **Neon ($Ne$)**: $x = 2.31$ 
- **Water ($H_{2}O$)**: $x = 0.52$ 
- **Benzene ($C_6H_6$)**: $x = 0.06$ 
- **Xenon ($Xe$)**: $x = 0.21$ (with small-core relativistic ECP)

*Note: The strict bounds enforced by this framework successfully detect and reject non-physical states (e.g., catching Hartree-Fock instabilities in Naphthalene when RHF finite-field perturbations fail to converge).*

### B. Future libMBD Comparison
**Next steps for validation**: We plan to directly compare the outputs of our explicitly bounded structural matrices against `libMBD`. A 1-to-1 numerical equivalence will firmly ground the formal Lean proofs into established production pipelines.

---

## Installation

You can globally install the numerical framework via PyPI:
```bash
pip install mbd-framework
```

## Global Command Line Usage

Once installed, the framework registers three native CLI endpoints.

### 1. Atomic Density Bounds Extraction (`mbd-compute`)
Computes the atomic polarizability arrays using PySCF and sets TS scaling parameters into `database.json`.
* **Example:**
  ```bash
  mbd-compute --molecule Benzene --basis aug-cc-pVDZ
  ```

### 2. Crystal Dispersion Simulation (`mbd-crystal`)
Resolves Cartesian lattice macroscopic dispersion scaling against empirical Pauli Repulsion logic.
* **Example:**
  ```bash
  mbd-crystal --target Benzene --epsilon 1.0
  ```

### 3. SERS Mathematical Equivalence (`mbd-sers`)
Tests numerical outputs comparing SERS exponential quenching strictly against MBD intrinsic boundaries.
* **Example:**
  ```bash
  mbd-sers --target Benzene --epsilon 2.0
  ```

## Repository Organization

- **`lean/`**
  - **`MBD_Theory.lean`**: Derivation of the $x = V_{\text{Bohr}} / \alpha$ proportionality limit.
  - **`Molecular_MBD.lean`**: Geometry matrix evaluations and `screened_is_bounded_unscreened` theorem.
- **`mbd_framework/`**
  - **`compute_volumes.py`**: PySCF derivation engine enforcing the Lean bounds.
  - **`crystal_validation.py`**: Atomic grid lattice calculator.
  - **`sers_unification.py`**: Analytical mapping logic for SERS parameters.
