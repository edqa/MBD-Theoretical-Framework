import numpy as np
import json
import os
import copy
from pyscf import gto, scf

def compute_molecule(name, geometry, basis='sto-3g'):
    print(f"\n--- Processing {name} ({basis}) ---")
    mol = gto.M(atom=geometry, basis=basis, verbose=3)
    
    # 1. SCF Energy
    mf = scf.RHF(mol).run()
    if not mf.converged:
        print(f"Warning: SCF did not converge for {name}")
    
    # 2. Polarizability Tensor (Finite-Field derivation)
    print("Computing Polarizability Tensor (Finite Field)...")
    h1 = mol.intor('int1e_r', comp=3)
    alpha_tensor = np.zeros((3,3))
    field = 1e-4
    hcore = mf.get_hcore()
    for i in range(3):
        mf_p = copy.copy(mf)
        mf_p.direct_scf = False
        mf_p.get_hcore = lambda *args, j=i, hc=hcore: hc + field * h1[j]
        mf_p.run(verbose=0)
        dip_p = mf_p.dip_moment(mol, mf_p.make_rdm1(), unit='A.U.', verbose=0)
        
        mf_m = copy.copy(mf)
        mf_m.direct_scf = False
        mf_m.get_hcore = lambda *args, j=i, hc=hcore: hc - field * h1[j]
        mf_m.run(verbose=0)
        dip_m = mf_m.dip_moment(mol, mf_m.make_rdm1(), unit='A.U.', verbose=0)
        
        alpha_tensor[i] = (dip_m - dip_p) / (2 * field)
        
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
    alpha_eff = max(1e-4, abs(alpha_iso))
    x = V_Bohr_au / alpha_eff
    
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
    
    results = {}
    for name, geo in geometries.items():
        if args.molecule != "all" and name.lower() != args.molecule.lower():
            continue
        
        # Safely enforce ECP specifically for massive core elements if generic aug basis requested
        b = 'def2-svp' if (name == 'Xe' and 'aug' in args.basis.lower()) else args.basis
        
        data = compute_molecule(name, geo, basis=b)
        results[name] = data

    db_path = os.path.join(os.path.dirname(__file__), "database.json")
    with open(db_path, "w") as f:
        json.dump(results, f, indent=4)
    print(f"\nSaved results to {db_path}")

if __name__ == "__main__":
    main()
