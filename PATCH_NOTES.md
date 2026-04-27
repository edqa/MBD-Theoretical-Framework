# Patch Notes — Sign-correction patch

## Summary

This patch fixes a sign error in the finite-field polarizability calculation
in `mbd_framework/compute_volumes.py` and a related issue with the Xe basis
substitution. The published `x` values are unchanged in magnitude; the
underlying `alpha_iso` values are now physically correct (positive).

## Changes

### 1. Sign of the finite-field polarizability (the actual bug)

In `compute_volumes.py`, the polarizability formula was inverted:

```python
# Before (wrong):
alpha_tensor[i] = (dip_m - dip_p) / (2 * field)

# After (correct):
alpha_tensor[i] = (dip_p - dip_m) / (2 * field)
```

Polarizability is defined as μᵢ = αᵢⱼ Eⱼ, so the induced dipole grows *with*
the field. The correct finite-difference formula is therefore
αᵢⱼ = [μᵢ(+E) − μᵢ(−E)] / (2E).

The old formula gave negative α for every species in the database, which is
unphysical for any closed-shell ground state. The negative values were masked
downstream by `alpha_eff = max(1e-4, abs(alpha_iso))`, which produced
correct-magnitude `x` values but published broken `alpha_iso` numbers.

### 2. Removed the `abs()` rescue

The line `alpha_eff = max(1e-4, abs(alpha_iso))` was a band-aid that hid the
sign bug above. With the sign fixed, it's no longer needed. It has been
replaced by an explicit `ValueError` if α ≤ 0, so future sign issues fail
loudly rather than silently producing junk in the database.

### 3. Xenon basis substitution

The original code silently substituted `def2-svp` for Xe whenever an `aug-`
basis was requested. `def2-svp` for Xe in PySCF is all-electron, non-
relativistic, and not really defensible for a 54-electron atom. The
substitution is now to `def2-tzvp` with the matching `def2-tzvp` ECP
explicitly set, and a `NOTE:` line is printed so the substitution is
reproducible from the log alone.

This is also a correction to the *physics*, not just the documentation:
α(Xe) shifts from 12.9 a.u. (broken, no ECP) to 20.3 a.u. (correct, with
small-core relativistic ECP), and `x` shifts from 0.32 to 0.21.

### 4. PySCF ≥ 2.10 compatibility

The original `copy.copy(mf)` + `mf_p.direct_scf = False` pattern fails on
PySCF ≥ 2.10 with `AttributeError: 'NoneType' object has no attribute 'get'`
because `_opt` is now lazily initialized. Replaced with a fresh `scf.RHF(mol)`
per finite-field step that uses the converged unperturbed density matrix as
its initial guess. This is faster *and* compatible with current PySCF.

### 5. Checkpointing and database merging

`compute_volumes.py` now writes `database.json` after every species rather
than only at the end of the full sweep, and merges with the existing database
rather than overwriting. This means partial runs accumulate, and you can run
species one at a time (`--molecule Naphthalene`) without losing prior results.

## Updated benchmark values

| Species   | Original α  | Patched α   | Original x | Patched x | Notes |
|-----------|-------------|-------------|------------|-----------|-------|
| He        | −1.2843     | **+1.2843** | 3.2615     | 3.2615    | sign only |
| Ne        | −1.8068     | **+1.8068** | 2.3184     | 2.3184    | sign only |
| H₂O       | −8.0609     | **+8.0609** | 0.5196     | 0.5196    | sign only |
| Benzene   | −68.5759    | **+68.5759**| 0.0611     | 0.0611    | sign only |
| Xe        | −12.903     | **+20.316** | 0.3246     | **0.2062**| sign + ECP |
| Naphthalene | (not run) | run locally | —          | —         | new |

α values are in atomic units (Bohr³). x = V_Bohr / α is dimensionless.

## Downstream files

`crystal_validation.py` and `sers_unification.py` consume only `data["x"]`
from the database, not `data["alpha_iso"]`, so the −131.99 kJ/mol benzene
lattice energy and the SERS equivalence check are unaffected by the sign
correction. Re-running them is recommended for completeness but no numerical
change is expected for He/Ne/H₂O/Benzene. The Xe row will change because of
the ECP substitution.

## Lean side

`MBDSystem` and `AtomNode` in the Lean files require `0 < alpha` as a
precondition. The original Python emitted negative α values, which means
those structures could not have been instantiated with the published numbers
— the Lean ↔ Python correspondence was technically broken at the type level.
With the sign fix, the database values now satisfy the Lean preconditions and
the formal-verification chain holds.

No changes to the Lean files are included in this patch. Two suggestions for
a future revision (lower priority):

- The headline theorem `screened_is_bounded_unscreened` reduces to the
  trivial property "for ε ≥ 1 and x > 0, ε^(−x) ≤ 1". Worth either
  strengthening to a non-trivial physical bound or framing more modestly in
  the README.
- `Molecular_MBD.lean` uses `0 < alpha_mol.det` as the well-formedness
  condition on the polarizability tensor. Positive determinant is weaker than
  positive-definiteness; a tensor with two negative and one positive
  eigenvalue still has positive determinant. Switch to `alpha_mol.PosDef` for
  a physically correct constraint.
