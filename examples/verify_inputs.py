#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick validation of molecule JSON input format"""
import json
import sys
import os
from pathlib import Path

os.environ['PYTHONIOENCODING'] = 'utf-8'

OK = "[OK]"
FAIL = "[FAIL]"
WARN = "[WARN]"


def validate_molecule_json(path: Path) -> bool:
    print(f"Validating: {path.name}")
    try:
        with open(path, encoding='utf-8') as f:
            mol = json.load(f)
    except Exception as e:
        print(f"  {FAIL} JSON parse error: {e}")
        return False

    required_keys = ["atoms"]
    for k in required_keys:
        if k not in mol:
            print(f"  {FAIL} Missing key: {k}")
            return False

    atoms = mol["atoms"]
    if not isinstance(atoms, list):
        print(f"  {FAIL} 'atoms' must be a list")
        return False

    n = len(atoms)
    pi_electrons = 0
    elements = set()
    for i, atom in enumerate(atoms):
        for k in ["element", "x", "y", "z"]:
            if k not in atom:
                print(f"  {FAIL} Atom {i} missing key: {k}")
                return False
        elements.add(atom["element"])
        if atom["element"] in ["C", "N", "O", "S", "P"]:
            pi_electrons += 1

    print(f"  {OK} Atoms: {n} (elements: {', '.join(sorted(elements))})")
    print(f"  {OK} pi-conjugated atoms: {pi_electrons} (pi-electrons: {pi_electrons})")
    print(f"  {OK} Hamiltonian matrix size: {n}x{n} = {n*n} elements")

    if "bonds" in mol:
        bonds = mol["bonds"]
        bond_ok = True
        for j, bond in enumerate(bonds):
            for k in ["atom1", "atom2"]:
                if k not in bond:
                    print(f"  {FAIL} Bond {j} missing key: {k}")
                    bond_ok = False
                    continue
                if not (0 <= bond[k] < n):
                    print(f"  {FAIL} Bond {j} {k}={bond[k]} out of range [0, {n-1}]")
                    bond_ok = False
        if bond_ok:
            print(f"  {OK} Bonds: {len(bonds)}")
    else:
        print(f"  {WARN} No bonds specified (will auto-infer from distances)")

    print(f"  {OK} Result: PASS\n")
    return True


def main():
    if sys.platform.startswith('win'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')

    project_root = Path(__file__).parent.parent
    mol_dir = project_root / "examples" / "molecules"

    print("=" * 60)
    print("  Huckel Engine - Input Validation")
    print("=" * 60 + "\n")

    all_ok = True
    json_files = sorted(mol_dir.glob("*.json"))
    if not json_files:
        print(f"{FAIL} No JSON molecule files found in {mol_dir}")
        sys.exit(1)

    for json_file in json_files:
        if not validate_molecule_json(json_file):
            all_ok = False

    print("=" * 60)
    if all_ok:
        print(f"All inputs valid! {OK}")
        print("\nNext steps:")
        print("  1. Install Rust: https://rustup.rs/")
        print("  2. Run: python run_huckel.py examples/molecules/benzene.json")
        print("  3. For GPU-accelerated 500+ atom systems:")
        print("     python run_huckel.py --generate-graphene 16x16 --grid-res 60")
    else:
        print(f"Some inputs have errors. {FAIL}")
        sys.exit(1)
    print("=" * 60)


if __name__ == "__main__":
    main()
