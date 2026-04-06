import json
import math
import os
import argparse
import numpy as np

def get_internal_coords(molecule="Benzene"):
    """Returns typical atomic (element, np.array([x,y,z])) for internal structure."""
    if molecule == "Benzene":
        return [
            ("C", np.array([0.00000,  1.40000,  0.00000])), ("C", np.array([1.21244,  0.70000,  0.00000])),
            ("C", np.array([1.21244, -0.70000,  0.00000])), ("C", np.array([0.00000, -1.40000,  0.00000])),
            ("C", np.array([-1.21244, -0.70000, 0.00000])), ("C", np.array([-1.21244,  0.70000, 0.00000])),
            ("H", np.array([0.00000,  2.49000,  0.00000])), ("H", np.array([2.15637,  1.24500,  0.00000])),
            ("H", np.array([2.15637, -1.24500,  0.00000])), ("H", np.array([0.00000, -2.49000,  0.00000])),
            ("H", np.array([-2.15637, -1.24500, 0.00000])), ("H", np.array([-2.15637,  1.24500, 0.00000]))
        ]
    elif molecule == "Naphthalene":
        return [
            ("C", np.array([-0.6976, 0.0, 1.2185])), ("C", np.array([0.6976, 0.0, 1.2185])),
            ("C", np.array([-1.3932, 0.0, 0.0])), ("C", np.array([1.3932, 0.0, 0.0])),
            ("C", np.array([-0.6976, 0.0, -1.2185])), ("C", np.array([0.6976, 0.0, -1.2185])),
            ("C", np.array([-2.7885, 0.0, 0.0])), ("C", np.array([2.7885, 0.0, 0.0])),
            ("C", np.array([-3.4862, 0.0, -1.2185])), ("C", np.array([3.4862, 0.0, -1.2185])),
            ("H", np.array([-1.2384, 0.0, 2.1549])), ("H", np.array([1.2384, 0.0, 2.1549])),
            ("H", np.array([-1.2384, 0.0, -2.1549])), ("H", np.array([1.2384, 0.0, -2.1549])),
            ("H", np.array([-3.3293, 0.0, 0.9364])), ("H", np.array([3.3293, 0.0, 0.9364])),
            ("H", np.array([-4.5714, 0.0, -1.2185])), ("H", np.array([4.5714, 0.0, -1.2185]))
        ]
    elif molecule == "Ice":
        return [
            ("O", np.array([0.0, 0.0, 0.0])),
            ("H", np.array([0.0, 0.757, 0.587])),
            ("H", np.array([0.0, -0.757, 0.587]))
        ]
    return []

def build_atomic_lattice(molecule="Benzene", a=7.46, b=9.67, c=7.03, supercell=3):
    """ Builds atomic-resolution molecular proxy crystal. """
    fractional_coords = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.0]),
        np.array([0.0, 0.5, 0.5]),
        np.array([0.5, 0.0, 0.5])
    ]
    internal_coords = get_internal_coords(molecule)
    
    lattice_molecules = []
    
    # Target origin molecule
    origin_mol = [(ele, pos) for ele, pos in internal_coords]
    
    for i in range(-supercell, supercell + 1):
        for j in range(-supercell, supercell + 1):
            for k in range(-supercell, supercell + 1):
                offset = np.array([i * a, j * b, k * c])
                
                for idx, fcoord in enumerate(fractional_coords):
                    # Skip the origin molecule to avoid internal interactions
                    if i == 0 and j == 0 and k == 0 and idx == 0:
                        continue
                        
                    mol_center = offset + fcoord * np.array([a, b, c])
                    molecule_atoms = [(ele, mol_center + pos) for ele, pos in internal_coords]
                    lattice_molecules.append(molecule_atoms)
                    
    return origin_mol, lattice_molecules

def compute_lattice_energy(origin_mol, lattice_molecules, x_mol, epsilon=1.0):
    bohr_to_ang = 0.529177
    hartree_to_kjmol = 2625.5
    
    # Generic unscaled empirical atom-atom C6 coefficients in atomic units (TS parameters)
    C6_dict = {
        ("C", "C"): 46.6,
        ("C", "H"): 17.5,
        ("H", "C"): 17.5,
        ("H", "H"): 6.5,
        ("O", "O"): 15.6,
        ("O", "H"): 10.0,
        ("H", "O"): 10.0
    }
    
    total_energy_au = 0.0
    screening_factor = math.pow(epsilon, -x_mol)
    
    # Extract structural arrays
    # Broadcast computation
    for ext_mol in lattice_molecules:
        for origin_ele, origin_pos in origin_mol:
            for ext_ele, ext_pos in ext_mol:
                
                dist_ang = np.linalg.norm(ext_pos - origin_pos)
                if dist_ang < 0.1: continue
                
                dist_bohr = dist_ang / bohr_to_ang
                c6_ij = C6_dict[(origin_ele, ext_ele)]
                
                v_ij = - (c6_ij * screening_factor) / (dist_bohr**6)
                total_energy_au += v_ij
                
    # Since we iterated 1 molecule against all others... The single mol interaction energy is accurate.
    # To get binding energy per molecule in crystal, div by 2
    e_lattice_kjmol = (total_energy_au / 2.0) * hartree_to_kjmol
    return e_lattice_kjmol

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon", type=float, default=1.0, help="Bulk Dielectric Constant")
    parser.add_argument("--target", type=str, default="Benzene", help="Target Molecule to compute (Benzene, Naphthalene, Ice)")
    args = parser.parse_args()

    db_path = os.path.join(os.getcwd(), "database.json")
    if not os.path.exists(db_path):
        print(f"Database not found at {db_path}! Run mbd-compute first.")
        return
        
    with open(db_path, "r") as f:
        db = json.load(f)
        
    target_db_key = "H2O" if args.target.lower() == "ice" else args.target
    molecule_key = next((k for k in db.keys() if k.lower() == target_db_key.lower()), None)
    if not molecule_key:
        print(f"{target_db_key} not computed in database.")
        return
        
    data = db[molecule_key]
    x_mol = data["x"]
    
    print(f"--- {args.target} Molecular Crystal Proxy ---")
    print(f"Computed intrinsic screening limit x = {x_mol:.4f}")
    
    print(f"\nBuilding 7x7x7 atomic proxy lattice for {args.target}...")
    origin_mol, lattice_molecules = build_atomic_lattice(molecule=args.target, supercell=3)
    print(f"Total Molecules in bounding box: {len(lattice_molecules) + 1}")
    print(f"Total Atoms evaluated iteratively: {len(origin_mol) * len(lattice_molecules)}")
    
    e_lattice = compute_lattice_energy(origin_mol, lattice_molecules, x_mol, epsilon=args.epsilon)
    
    print(f"\n--- Lattice Energy Calculation ---")
    print(f"Dielectric Constant (ε_eff) = {args.epsilon}")
    print(f"MBD Screening Factor (ε^-x) = {math.pow(args.epsilon, -x_mol):.5f}")
    print(f"Computed Dispersion Energy = {e_lattice:.2f} kJ/mol")

if __name__ == "__main__":
    main()
