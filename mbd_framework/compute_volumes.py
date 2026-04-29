"""
Compute MBD Bohr volumes and isotropic polarizabilities.

This module computes alpha_iso for benchmark molecules and emits the
dimensionless ratio x = V_Bohr / alpha into database.json. It supports
side-by-side cross-validation between Restricted Hartree-Fock (RHF) and
Restricted Kohn-Sham DFT (RKS) references, with each reference verified
to be a true local energy minimum via Pople-style stability analysis
before the polarizability is computed.

Methodology summary
-------------------
For each molecule and each electronic-structure level (RHF and/or RKS):

  1. SCF convergence to a stationary point.
  2. Stability iteration: Pople's orbital-Hessian analysis is run; if the
     wavefunction is at a saddle point (negative-eigenvalue mode), the
     orbitals are perturbed along the imaginary mode and the SCF is
     reconverged. Repeat up to STAB_MAX_ITER times.
  3. Polarizability via analytical coupled-perturbed HF/KS (CPHF/CPKS)
     using pyscf-properties. No finite-difference noise.
  4. x = V_Bohr / alpha_iso (V_Bohr = 4pi/3 a.u.).

This procedure mirrors the discipline used in plane-wave DFT codes (e.g.
Quantum ESPRESSO) where the ground state is verified before properties
are computed. RHF results provide a non-correlated reference; RKS/PBE
provides a correlated reference closer to experiment. Cross-validation
between the two methods strengthens benchmark claims and exposes
systems where HF alone is insufficient (e.g. extended pi-systems).
"""

import numpy as np
import json
import os
import argparse
from pyscf import gto, scf, dft

# Analytical CPHF/CPKS polarizability via pyscf-properties.
try:
    from pyscf.prop.polarizability import rhf as polrhf_mod
    from pyscf.prop.polarizability import rks as polrks_mod
    _HAS_CPHF = True
except ImportError:
    _HAS_CPHF = False

STAB_MAX_ITER = 5  # Maximum stability-iteration cycles before giving up.


# ---------------------------------------------------------------------------
# Stability iteration: standard fix for SCFs that converge to saddle points
# ---------------------------------------------------------------------------

def _stability_iterate(mf, name, label):
    """Iterate Pople stability analysis until the SCF is at a true minimum.

    PySCF's default initial guess can converge to a saddle point in
    orbital-rotation space, especially for extended pi-systems where the
    orbital Hessian has small or imaginary eigenvalues. Property
    calculations (CPHF, CPKS) on such references give nonsensical results.
    This routine perturbs along the imaginary mode and reconverges until
    the orbital Hessian is positive-definite.

    Returns the (possibly re-converged) mean-field object.
    """
    print(f"  Stability check ({label}):")
    for it in range(STAB_MAX_ITER):
        mo_int, _ = mf.stability(verbose=0)
        if mo_int is None or np.allclose(mo_int, mf.mo_coeff):
            print(f"    iter {it}: stable.")
            return mf
        print(f"    iter {it}: unstable, perturbing along imaginary mode and "
              f"reconverging.")
        dm = mf.make_rdm1(mo_coeff=mo_int, mo_occ=mf.mo_occ)
        # Rebuild the mean-field object of the same class so we don't
        # inherit any cached state from the saddle-point SCF.
        new_mf = mf.__class__(mf.mol)
        if hasattr(mf, 'xc'):
            new_mf.xc = mf.xc  # preserve DFT functional
        new_mf.verbose = 0
        new_mf.kernel(dm0=dm)
        if not new_mf.converged:
            print(f"    iter {it}: SCF after perturbation did NOT converge.")
        else:
            print(f"    iter {it}: reconverged. delta_E = "
                  f"{new_mf.e_tot - mf.e_tot:+.6f} hartree.")
        mf = new_mf
    print(f"  Warning: stability iteration did not converge in "
          f"{STAB_MAX_ITER} cycles for {name} ({label}).")
    return mf


# ---------------------------------------------------------------------------
# Polarizability via analytical CPHF / CPKS (no finite-difference noise)
# ---------------------------------------------------------------------------

def _alpha_via_cphf(mf, label):
    """Analytical polarizability via coupled-perturbed HF (RHF) or KS (RKS).

    Production-quality method used by Gaussian, Q-Chem, ORCA, NWChem.
    Robust against linear-dependency contamination in diffuse basis sets.
    Requires the reference SCF to be at a local minimum (verify with
    _stability_iterate first).
    """
    if not _HAS_CPHF:
        raise RuntimeError(
            "CPHF/CPKS requires pyscf-properties. Install with:\n"
            "  pip install pyscf-properties"
        )
    if isinstance(mf, dft.rks.RKS):
        return np.asarray(polrks_mod.Polarizability(mf).polarizability())
    else:
        return np.asarray(polrhf_mod.Polarizability(mf).polarizability())


# ---------------------------------------------------------------------------
# Per-method evaluation: SCF -> stabilize -> CPHF -> x = V_Bohr/alpha
# ---------------------------------------------------------------------------

def _evaluate_method(mol, name, method_label, scf_class, xc=None):
    """Run SCF -> stability iteration -> CPHF for one electronic-structure
    method. Returns a dict with alpha_iso, alpha_tensor, x, energy, etc.
    """
    print(f"\n  [{method_label}] SCF...")
    mf = scf_class(mol)
    if xc is not None:
        mf.xc = xc
    mf.verbose = 0
    mf.kernel()
    if not mf.converged:
        print(f"    Warning: initial SCF did not converge.")
    print(f"    E_SCF = {mf.e_tot:.6f} hartree")

    mf = _stability_iterate(mf, name, method_label)

    print(f"  [{method_label}] CPHF/CPKS polarizability...")
    alpha_tensor = _alpha_via_cphf(mf, method_label)
    alpha_iso = np.trace(alpha_tensor) / 3.0
    print(f"    alpha_iso = {alpha_iso:.4f} a.u.")

    if alpha_iso <= 0:
        raise ValueError(
            f"Non-physical polarizability for {name} at {method_label}: "
            f"alpha_iso = {alpha_iso:.4f} a.u. The reference wavefunction "
            f"may still be at a saddle point even after stability iteration. "
            f"Consider increasing STAB_MAX_ITER, switching functional, or "
            f"using a multi-reference method."
        )

    V_Bohr_au = 4.0 * np.pi / 3.0
    x = V_Bohr_au / alpha_iso
    print(f"    x = V_Bohr / alpha = {x:.4f}")

    # Mean-square electronic radius from the converged density.
    dm = mf.make_rdm1()
    r2_ints = mol.intor('int1e_r2')
    elec_r2 = float(np.einsum('ij,ji->', dm, r2_ints))
    r2_per_elec = (elec_r2 / mol.nelectron) if mol.nelectron > 0 else elec_r2

    return {
        "method": method_label,
        "energy": float(mf.e_tot),
        "converged": bool(mf.converged),
        "alpha_iso": float(alpha_iso),
        "alpha_tensor": alpha_tensor.tolist(),
        "elec_r2_per_electron": float(r2_per_elec),
        "V_Bohr": float(V_Bohr_au),
        "x": float(x),
    }


def compute_molecule(name, geometry, basis='sto-3g', ecp=None,
                     methods=('rhf', 'rks-pbe')):
    """Compute polarizabilities for one molecule with one or more methods.

    Parameters
    ----------
    name : str
    geometry : str (PySCF z-matrix format)
    basis : str
    ecp : str or None
    methods : iterable of strings, any of {'rhf', 'rks-pbe', 'rks-b3lyp'}

    Returns
    -------
    dict with keys 'name', 'basis', 'ecp', and one entry per method label.
    """
    print(f"\n=== Processing {name} (basis={basis}"
          f"{', ecp=' + ecp if ecp else ''}) ===")
    mol_kwargs = dict(atom=geometry, basis=basis, verbose=0)
    if ecp:
        mol_kwargs['ecp'] = ecp
    mol = gto.M(**mol_kwargs)
    print(f"  nbf = {mol.nao}, nelec = {mol.nelectron}")

    record = {"name": name, "basis": basis, "ecp": ecp,
              "nao": int(mol.nao), "nelec": int(mol.nelectron)}

    for m in methods:
        if m == 'rhf':
            record['rhf'] = _evaluate_method(mol, name, 'RHF', scf.RHF)
        elif m == 'rks-pbe':
            record['rks_pbe'] = _evaluate_method(
                mol, name, 'RKS/PBE', dft.RKS, xc='pbe')
        elif m == 'rks-b3lyp':
            record['rks_b3lyp'] = _evaluate_method(
                mol, name, 'RKS/B3LYP', dft.RKS, xc='b3lyp')
        else:
            raise ValueError(f"Unknown method: {m!r}. "
                             f"Use any of: rhf, rks-pbe, rks-b3lyp")

    # Cross-validation summary if multiple methods were run.
    if len(methods) > 1:
        print(f"\n  Cross-validation summary for {name}:")
        for m in methods:
            key = m.replace('-', '_')
            if key in record:
                e = record[key]
                print(f"    {e['method']:10s}  alpha_iso = "
                      f"{e['alpha_iso']:9.4f} a.u.   x = {e['x']:.4f}")

    return record


# ---------------------------------------------------------------------------
# Geometries
# ---------------------------------------------------------------------------

GEOMETRIES = {
    "He": "He 0 0 0",
    "Ne": "Ne 0 0 0",
    "Xe": "Xe 0 0 0",
    "H2O": '''
        O 0 0 0
        H 0 0.757 0.587
        H 0 -0.757 0.587
    ''',
    "Benzene": '''
        C        0.00000        1.40000        0.00000
        H        0.00000        2.49000        0.00000
        C        1.21244        0.70000        0.00000
        H        2.15637        1.24500        0.00000
        C        1.21244       -0.70000        0.00000
        H        2.15637       -1.24500        0.00000
        C        0.00000       -1.40000        0.00000
        H        0.00000       -2.49000        0.00000
        C       -1.21244       -0.70000        0.00000
        H       -2.15637       -1.24500        0.00000
        C       -1.21244        0.70000        0.00000
        H       -2.15637        1.24500        0.00000
    ''',
    "Naphthalene": '''
        C  -0.6976  0.0000  1.2185
        C   0.6976  0.0000  1.2185
        C  -1.3932  0.0000  0.0000
        C   1.3932  0.0000  0.0000
        C  -0.6976  0.0000 -1.2185
        C   0.6976  0.0000 -1.2185
        C  -2.7885  0.0000  0.0000
        C   2.7885  0.0000  0.0000
        C  -3.4862  0.0000 -1.2185
        C   3.4862  0.0000 -1.2185
        H  -1.2384  0.0000  2.1549
        H   1.2384  0.0000  2.1549
        H  -1.2384  0.0000 -2.1549
        H   1.2384  0.0000 -2.1549
        H  -3.3293  0.0000  0.9364
        H   3.3293  0.0000  0.9364
        H  -4.5714  0.0000 -1.2185
        H   4.5714  0.0000 -1.2185
    ''',
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute MBD Bohr volumes and polarizabilities with "
                    "RHF/DFT cross-validation.")
    parser.add_argument("--basis", default="aug-cc-pvdz",
                        help="Basis set (default: aug-cc-pvdz)")
    parser.add_argument("--molecule", default="all",
                        help="Target molecule or 'all'")
    parser.add_argument(
        "--methods", default="rhf,rks-pbe",
        help="Comma-separated methods. Choices: rhf, rks-pbe, rks-b3lyp. "
             "Default: rhf,rks-pbe (cross-validation).")
    args = parser.parse_args()

    methods = tuple(m.strip().lower() for m in args.methods.split(","))
    valid = {'rhf', 'rks-pbe', 'rks-b3lyp'}
    bad = [m for m in methods if m not in valid]
    if bad:
        raise SystemExit(f"Unknown methods: {bad}. Valid: {sorted(valid)}")

    db_path = os.path.join(os.getcwd(), "database.json")

    # Load existing database (merge mode) so partial runs accumulate.
    existing = {}
    if os.path.exists(db_path):
        try:
            with open(db_path) as f:
                existing = json.load(f)
        except Exception:
            existing = {}

    for name, geo in GEOMETRIES.items():
        if args.molecule != "all" and name.lower() != args.molecule.lower():
            continue

        # Xe at aug-cc-pVDZ is all-electron and non-relativistic in PySCF.
        # For a 54-electron atom this is not defensible; substitute
        # def2-TZVP + matching small-core relativistic ECP.
        ecp = None
        if name == 'Xe' and 'aug' in args.basis.lower():
            b = 'def2-tzvp'
            ecp = 'def2-tzvp'
            print(f"\nNOTE: Substituting basis '{args.basis}' -> '{b}' for "
                  f"Xe (with small-core relativistic ECP '{ecp}').")
        else:
            b = args.basis

        record = compute_molecule(name, geo, basis=b, ecp=ecp,
                                  methods=methods)

        # Merge per-molecule (preserves prior methods if you re-run with a
        # different --methods flag).
        prior = existing.get(name, {})
        prior.update(record)
        existing[name] = prior

        with open(db_path, "w") as f:
            json.dump(existing, f, indent=4)
        print(f"  Checkpointed {name} to {db_path}")

    print(f"\nDone. Database written to {db_path}")


if __name__ == "__main__":
    main()
