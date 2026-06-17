import json
import math
import sys
import os

def generate_graphene_sheet(n_x: int, n_y: int, filename: str):
    a_cc = 1.42
    a1_x = math.sqrt(3) * a_cc
    a2_y = 3.0 * a_cc

    atoms = []
    bonds = []
    atom_idx = 0

    bond_set = set()

    def add_bond(i, j):
        if i > j:
            i, j = j, i
        key = (i, j)
        if key not in bond_set:
            bond_set.add(key)
            bonds.append({"atom1": i, "atom2": j, "bond_order": 1.33})

    atom_map = {}

    for i_y in range(n_y):
        for i_x in range(n_x):
            for sub in range(2):
                x = i_x * a1_x + (sub % 2) * (a1_x / 2)
                y = i_y * a2_y + (sub % 2) * (a_cc * 1.5)
                z = 0.0

                atoms.append({
                    "element": "C",
                    "x": round(x, 6),
                    "y": round(y, 6),
                    "z": round(z, 6)
                })
                atom_map[(i_y, i_x, sub)] = atom_idx
                atom_idx += 1

    for i_y in range(n_y):
        for i_x in range(n_x):
            for sub in range(2):
                idx = atom_map.get((i_y, i_x, sub))
                if idx is None:
                    continue

                if sub == 0:
                    right = atom_map.get((i_y, i_x, 1))
                    if right is not None:
                        add_bond(idx, right)

                    right_up = atom_map.get((i_y + 1, i_x, 1))
                    if right_up is not None:
                        add_bond(idx, right_up)

                    left_up = atom_map.get((i_y + 1, i_x - 1, 1))
                    if left_up is not None:
                        add_bond(idx, left_up)

                else:
                    left = atom_map.get((i_y, i_x, 0))
                    if left is not None:
                        add_bond(idx, left)

                    right_down = atom_map.get((i_y - 1, i_x, 0))
                    if right_down is not None:
                        add_bond(idx, right_down)

                    left_down = atom_map.get((i_y - 1, i_x + 1, 0))
                    if left_down is not None:
                        add_bond(idx, left_down)

    h_atom_start = atom_idx
    bond_set_copy = set()
    for (i, j) in bond_set:
        bond_set_copy.add((i, j))
        bond_set_copy.add((j, i))

    for i_y in range(n_y):
        for i_x in range(n_x):
            for sub in range(2):
                idx = atom_map.get((i_y, i_x, sub))
                if idx is None:
                    continue

                neighbor_count = 0
                for (a, b) in bond_set:
                    if a == idx or b == idx:
                        neighbor_count += 1

                if neighbor_count < 3:
                    atom = atoms[idx]
                    dx = 0.0
                    dy = 0.0
                    neighbors = []
                    for (a, b) in bond_set:
                        if a == idx:
                            neighbors.append(b)
                        elif b == idx:
                            neighbors.append(a)

                    for n_idx in neighbors:
                        dx += atoms[n_idx]["x"] - atom["x"]
                        dy += atoms[n_idx]["y"] - atom["y"]

                    norm = math.sqrt(dx * dx + dy * dy)
                    if norm > 1e-10:
                        dx /= norm
                        dy /= norm
                    else:
                        dx = -1.0 if sub == 1 else 1.0
                        dy = 0.0

                    h_x = atom["x"] - dx * 1.09
                    h_y = atom["y"] - dy * 1.09
                    h_z = 0.0

                    atoms.append({
                        "element": "H",
                        "x": round(h_x, 6),
                        "y": round(h_y, 6),
                        "z": round(h_z, 6)
                    })
                    bonds.append({
                        "atom1": idx,
                        "atom2": h_atom_start,
                        "bond_order": 1.0
                    })
                    h_atom_start += 1

    molecule = {
        "name": f"Graphene Sheet {n_x}x{n_y} (C{2*n_x*n_y})",
        "charge": 0,
        "multiplicity": 1,
        "atoms": atoms,
        "bonds": bonds
    }

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as f:
        json.dump(molecule, f, indent=2)

    print(f"Generated: {filename}")
    print(f"  Atoms: {len(atoms)}")
    print(f"  Bonds: {len(bonds)}")
    print(f"  π-conjugated atoms (C): {2*n_x*n_y}")
    print(f"  Hamiltonian size: {2*n_x*n_y}x{2*n_x*n_y}")


if __name__ == "__main__":
    nx = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    ny = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(__file__), f"graphene_{nx}x{ny}.json")
    generate_graphene_sheet(nx, ny, out)
