import json
import os
import math

def main():
    db_path = os.path.join(os.path.dirname(__file__), "database.json")
    if not os.path.exists(db_path):
        print(f"Error: {db_path} not found.")
        return
        
    with open(db_path, "r") as f:
        database = json.load(f)
        
    print(f"Loaded {len(database)} systems from database.")
    
    # Test cases for bulk dielectric
    epsilons = [1.0, 2.0, 10.0, 80.0] # Vacuum, non-polar liquid, polar, water
    
    pairs = [
        ("He", "He"),
        ("Ne", "Ne"),
        ("Xe", "Xe"),
        ("H2O", "H2O"),
        ("Benzene", "Benzene"),
        ("Benzene", "H2O")
    ]
    
    for mol1, mol2 in pairs:
        if mol1 not in database or mol2 not in database:
            continue
            
        d1 = database[mol1]
        d2 = database[mol2]
        
        # Effective exponent x_ij = (x_i + x_j) / 2
        x1 = d1["x"]
        x2 = d2["x"]
        x_ij = (x1 + x2) / 2.0
        
        print(f"\n--- Interaction: {mol1} -- {mol2} ---")
        print(f"[{mol1}] x = {x1:.4f}")
        print(f"[{mol2}] x = {x2:.4f}")
        print(f"Average path screening exponent x_ij = {x_ij:.4f}")
        
        for eps in epsilons:
            # Pairwise geometric mean path screening factor = eps^(-x_ij)
            s_factor = math.pow(eps, -x_ij)
            print(f"  Dielectric ε={eps:<4}  =>  Screening Factor (ε^{{-x}}) = {s_factor:.5f}")
            
if __name__ == "__main__":
    main()
