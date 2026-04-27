import numpy as np
import json
import os
import copy
from pyscf import gto, scf

def compute_molecule(name, geometry, basis='sto-3g', ecp=None):
    print(f"\n--- Processing {name} ({basis}{', ecp=' + ecp if ecp else ''}) ---")
    mol_kwargs = dict(atom=geometry, basis=basis, verbose=3)
    if ecp:
        mol_kwargs['ecp'] = ecp
    mol = gto.M(**mol_kwargs)
    
    # 1. SCF Energy
    mf = scf.RHF(mol).run()
    if not mf.converged:
        print(f"Warning: SCF did not converge for {name}")
    
    # 2. Polarizability Tensor (Finite-Field derivation)
    # We rebuild a fresh RHF object for each ±field calculation rather than
    # mutating the converged one. The copy.copy/direct_scf=False pattern fails
    # on PySCF >= 2.10 because of changes to how _opt is lazily initialized.
    print("Computing Polarizability Tensor (Finite Field)...")
    h1 = mol.intor('int1e_r', comp=3)
    hcore0 = scf.hf.get_hcore(mol)
    alpha_tensor = np.zeros((3, 3))
    field = 1e-4

    def run_with_field(field_vec, dm0):
        mf_f = scf.RHF(mol)
        # Override get_hcore to add the dipole perturbation -E . r
        # (sign convention: H' = -mu . E = +e r . E for electrons)
        mf_f.get_hcore = lambda *args, **kwargs: (
            hcore0 + field_vec[0] * h1[0] + field_vec[1] * h1[1] + field_vec[2] * h1[2]
        )
        mf_f.verbose = 0
        # Start from converged unperturbed DM — much faster than atomic guess
        mf_f.kernel(dm0=dm0)
        return mf_f.dip_moment(mol, mf_f.make_rdm1(), unit='A.U.', verbose=0)

    dm_ref = mf.make_rdm1()
    for i in range(3):
        e_plus = np.zeros(3); e_plus[i] = +field
        e_minus = np.zeros(3); e_minus[i] = -field
        dip_p = run_with_field(e_plus, dm_ref)
        dip_m = run_with_field(e_minus, dm_ref)
        # Polarizability: mu_i = alpha_ij * E_j, so induced dipole grows WITH field.
        # Therefore alpha_ij = (mu(+E) - mu(-E)) / (2*E), NOT (mu(-E) - mu(+E)).
        alpha_tensor[i] = (dip_p - dip_m) / (2 * field)
        
    alpha_iso = np.trace(alpha_tensor) / 3.0
    print(f"Isotropic Polarizability: {alpha_iso:.4f} a.u.")
    
    # The mean-square radius of the total electron density 
    dm = mf.make_rdm1()
    r2_ints = mol.intor('int1e_r2')
    elec_r2 = float(np.einsum('ij,ji->', dm, r2_ints))
    r2_per_elec = elec_r2 / mol.nelectron if mol.nelectron > 0 else elec_r2
    
    # Fundamental Bohr Volume (Universal Constant: 4/3 * pi * a_0^3)
    # Using 1 bohr radius universally to map polarizability across the vacuum limit
    V_Bohr_au = (4.0 * np.pi / 3.0) 
    
    print(f"Universal V_Bohr: {V_Bohr_au:.4f} a.u. (~0.62 A^3)")
    
    # Theoretical x = V_Bohr / alpha
    # With the sign convention fixed, alpha_iso must be positive for any
    # well-behaved closed-shell ground state. Fail loudly if not.
    if alpha_iso <= 0:
        raise ValueError(
            f"Non-physical polarizability for {name}: alpha_iso = {alpha_iso:.4f} a.u. "
            f"Expected alpha_iso > 0 for a closed-shell ground state. "
            f"Check basis set, SCF convergence, and finite-field sign convention."
        )
    x = V_Bohr_au / alpha_iso
    
    print(f"Computed x = V_Bohr / alpha = {x:.4f}")
    
    return {
        "name": name,
        "alpha_iso": float(alpha_iso),
        "alpha_tensor": alpha_tensor.tolist(),
        "elec_r2": float(r2_per_elec),
        "V_Bohr": float(V_Bohr_au),
        "x": float(x)
    }

import argparse

def main():
    parser = argparse.ArgumentParser(description="Compute MBD Bohr volumes and polarizabilities.")
    parser.add_argument("--basis", type=str, default="aug-cc-pvdz", help="Basis set (e.g. sto-3g, aug-cc-pVTZ)")
    parser.add_argument("--molecule", type=str, default="all", help="Target molecule to calculate or 'all'")
    args = parser.parse_args()

    geometries = {
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
        '''
    }
    
    db_path = os.path.join(os.getcwd(), "database.json")

    def save_checkpoint():
        existing = {}
        if os.path.exists(db_path):
            try:
                with open(db_path) as f:
                    existing = json.load(f)
            except Exception:
                existing = {}
        existing.update(results)
        with open(db_path, "w") as f:
            json.dump(existing, f, indent=4)

    results = {}
    for name, geo in geometries.items():
        if args.molecule != "all" and name.lower() != args.molecule.lower():
            continue

        # Xe is a 54-electron atom and aug-cc-pVDZ for Xe in PySCF is not the
        # right reference — it's all-electron and lacks scalar-relativistic
        # corrections. Use def2-TZVP + matching small-core relativistic ECP.
        # This substitution is logged so the output is reproducible.
        ecp = None
        if name == 'Xe' and 'aug' in args.basis.lower():
            b = 'def2-tzvp'
            ecp = 'def2-tzvp'
            print(f"NOTE: Substituting basis '{args.basis}' -> '{b}' for Xe "
                  f"(with small-core relativistic ECP '{ecp}').")
        else:
            b = args.basis

        data = compute_molecule(name, geo, basis=b, ecp=ecp)
        results[name] = data
        save_checkpoint()
        print(f"Checkpointed {name} to {db_path}")

    print(f"\nSaved results to {db_path}")

if __name__ == "__main__":
    main()
