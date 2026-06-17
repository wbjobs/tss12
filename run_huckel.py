#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hückel Quantum Mechanics Calculation Runner

通过PyO3调用Rust库执行半经验量子力学计算。

用法:
    python run_huckel.py <molecule.json> [--output DIR] [--no-gpu]
                         [--grid-res N] [--grid-padding Å]
                         [--generate-graphene NxN]
                         [--force-algorithm simple|ehmo|auto]

示例:
    python run_huckel.py examples/molecules/benzene.json
    python run_huckel.py examples/molecules/naphthalene.json --output results/naphthalene
    python run_huckel.py --generate-graphene 16x16 --grid-res 64
    python run_huckel.py examples/molecules/benzene.json --no-gpu
    python run_huckel.py examples/molecules/ptcl4.json --force-algorithm ehmo
"""

import argparse
import json
import os
import sys
import subprocess
import importlib
import importlib.util
import time
import webbrowser
from pathlib import Path


ATOMIC_NUMBER_TABLE = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92,
}

TRANSITION_METAL_RANGES = [
    (21, 30), (39, 48), (57, 80), (89, 103),
]

LANTHANIDE_RANGE = (57, 71)
ACTINIDE_RANGE = (89, 103)


def get_atomic_number(element_symbol: str) -> int:
    """Return atomic number for element symbol (case-insensitive)."""
    if not element_symbol:
        return 0
    key = element_symbol[0].upper() + element_symbol[1:].lower()
    return ATOMIC_NUMBER_TABLE.get(key, 0)


def is_heavy_metal(element_symbol: str) -> bool:
    """Detect if an element is a heavy/transition metal requiring EHMO."""
    z = get_atomic_number(element_symbol)
    if z == 0:
        return False
    for start, end in TRANSITION_METAL_RANGES:
        if start <= z <= end:
            return True
    if LANTHANIDE_RANGE[0] <= z <= LANTHANIDE_RANGE[1]:
        return True
    if ACTINIDE_RANGE[0] <= z <= ACTINIDE_RANGE[1]:
        return True
    if z >= 78:
        return True
    return False


def contains_heavy_metal(mol_data: dict) -> list:
    """Scan atoms list for heavy metals, return list of (element, Z)."""
    heavy = []
    atoms = mol_data.get("atoms", [])
    seen = set()
    for atom in atoms:
        elem = atom.get("element", "")
        if elem in seen:
            continue
        seen.add(elem)
        z = get_atomic_number(elem)
        if is_heavy_metal(elem):
            heavy.append((elem, z))
    return heavy


def detect_algorithm(mol_data: dict, user_force: str = "auto") -> tuple:
    """
    Auto-detect which algorithm to use.
    Returns (algorithm_name, reason_english, heavy_metals_list)
    """
    heavy = contains_heavy_metal(mol_data)

    if user_force and user_force.lower() in ("simple", "simplehuckel", "huckel"):
        if heavy:
            reason = f"[WARN] User forced Simple Hückel, but heavy metals detected: {[h[0] for h in heavy]}! Consider --force-algorithm ehmo"
        else:
            reason = "Algorithm set to Simple Hückel by user."
        return "simple", reason, heavy

    if user_force and user_force.lower() in ("ehmo", "extended", "extendedhuckel"):
        return "ehmo", "Algorithm set to Extended Hückel (EHMO) by user.", heavy

    if heavy:
        names = ", ".join(f"{e}(Z={z})" for e, z in heavy)
        return "ehmo", f"Heavy metals detected [{names}] -> Extended Hückel (EHMO) selected.", heavy

    return "simple", "Only light elements (H/C/N/O/S/P) -> Simple Hückel selected.", heavy


def check_rust_installed() -> bool:
    """检查Rust工具链是否可用"""
    try:
        result = subprocess.run(
            ["rustc", "--version"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            print(f"[Setup] Rust detected: {result.stdout.strip()}")
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    try:
        subprocess.run(
            ["cargo", "--version"],
            capture_output=True, text=True, timeout=10
        )
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def find_cargo() -> str:
    """尝试在常见位置查找cargo"""
    cargo_paths = [
        os.path.expanduser("~/.cargo/bin/cargo"),
        os.path.expanduser("~/.cargo/bin/cargo.exe"),
        "cargo",
    ]
    for p in cargo_paths:
        try:
            r = subprocess.run([p, "--version"], capture_output=True, text=True, timeout=5)
            if r.returncode == 0:
                return p
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return None


def build_extension(project_root: Path) -> bool:
    """使用maturin或setuptools-rust编译Python扩展"""
    print("[Setup] Building Rust extension module...")
    cargo = find_cargo()
    if not cargo:
        print("[Error] Cargo (Rust build tool) not found.")
        print("        Please install Rust via: https://rustup.rs/")
        return False

    try:
        import maturin  # type: ignore
        print("[Setup] Using maturin to build...")
        result = subprocess.run(
            [sys.executable, "-m", "maturin", "develop", "--release"],
            cwd=project_root,
            capture_output=False
        )
        return result.returncode == 0
    except ImportError:
        pass

    print("[Setup] Using cargo build + manual setup...")
    result = subprocess.run(
        [cargo, "build", "--release"],
        cwd=project_root,
        capture_output=False
    )
    if result.returncode != 0:
        print("[Error] Cargo build failed.")
        return False

    lib_ext = None
    lib_name = None
    if sys.platform.startswith("win"):
        lib_ext = "dll"
        lib_name = "huckel_engine.dll"
    elif sys.platform == "darwin":
        lib_ext = "dylib"
        lib_name = "libhuckel_engine.dylib"
    else:
        lib_ext = "so"
        lib_name = "libhuckel_engine.so"

    src_lib = project_root / "target" / "release" / lib_name
    if not src_lib.exists():
        print(f"[Error] Built library not found at {src_lib}")
        return False

    import shutil

    ext_suffix = ""
    try:
        import sysconfig
        ext_suffix = sysconfig.get_config_var("EXT_SUFFIX") or f".{lib_ext}"
    except Exception:
        ext_suffix = f".{lib_ext}"

    target_lib = project_root / f"huckel_engine{ext_suffix}"
    shutil.copy2(src_lib, target_lib)
    print(f"[Setup] Installed extension: {target_lib}")
    return True


def load_module(project_root: Path):
    """尝试加载Python扩展模块"""
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    try:
        import huckel_engine  # type: ignore
        print("[Setup] Extension module loaded successfully.")
        return huckel_engine
    except ImportError as e:
        print(f"[Setup] Module not yet built: {e}")
        return None


def generate_graphene(nx: int, ny: int, project_root: Path) -> Path:
    """生成石墨烯纳米片JSON文件"""
    mol_dir = project_root / "examples" / "molecules"
    mol_dir.mkdir(parents=True, exist_ok=True)
    gen_script = mol_dir / "generate_graphene.py"

    result_file = mol_dir / f"graphene_{nx}x{ny}.json"
    subprocess.run(
        [sys.executable, str(gen_script), str(nx), str(ny), str(result_file)],
        check=True
    )
    return result_file


def run_reaction_path_mode(args, module, project_root):
    """反应路径模拟模式"""
    reaction_path = Path(args.reaction)
    if not reaction_path.is_absolute():
        reaction_path = (project_root / reaction_path).resolve()

    if not reaction_path.exists():
        print(f"[Error] Reaction JSON not found: {reaction_path}")
        sys.exit(1)

    print(f"\n[Input] Reaction file: {reaction_path}")

    try:
        with open(reaction_path, "r", encoding="utf-8") as f:
            reaction_data = json.load(f)
    except Exception as e:
        print(f"[Error] Cannot parse reaction JSON: {e}")
        sys.exit(1)

    reactant = reaction_data.get("reactant", {})
    product = reaction_data.get("product", {})

    r_atoms = reactant.get("atoms", [])
    p_atoms = product.get("atoms", [])
    print(f"  Reactant: {len(r_atoms)} atoms ({', '.join(set(a['element'] for a in r_atoms))})")
    print(f"  Product:  {len(p_atoms)} atoms ({', '.join(set(a['element'] for a in p_atoms))})")
    print(f"  Steps:    {args.num_steps}")

    all_atoms = r_atoms + p_atoms
    _, reason, heavy = detect_algorithm({"atoms": all_atoms}, args.force_algorithm)
    print(f"  Algorithm: {reason}")
    if heavy:
        print(f"  Heavy metals: {', '.join(f'{e}(Z={z})' for e,z in heavy)}")

    output_dir = Path(args.output) if args.output else project_root / "results" / reaction_path.stem
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[Output] Results directory: {output_dir}")
    print(f"[Params] GPU: {not args.no_gpu}, Grid res: {args.grid_res}, Steps: {args.num_steps}")

    t0 = time.time()
    try:
        reaction_json_text = json.dumps(reaction_data)
        engine_force = None
        if args.force_algorithm == "ehmo":
            engine_force = "ehmo"
        elif args.force_algorithm == "simple":
            engine_force = "simple"

        engine = module.HuckelEngine(use_gpu=not args.no_gpu, force_algorithm=engine_force)
        result = engine.calculate_reaction_path(
            reaction_json_text,
            str(output_dir),
            num_steps=args.num_steps,
            grid_resolution=args.grid_res,
            grid_padding=args.grid_padding
        )
    except Exception as e:
        print(f"\n[Error] Reaction path simulation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"\n[Done] Total elapsed: {elapsed:.2f} seconds")

    html_path = output_dir / "reaction_path.html"
    json_path = output_dir / "reaction_path.json"
    for f in [json_path, html_path]:
        if f.exists():
            size_kb = f.stat().st_size / 1024
            print(f"  [OK] {f.relative_to(project_root)}  ({size_kb:.1f} KB)")

    if result:
        energies = result.get("total_energies", [])
        gaps = result.get("homo_lumo_gaps", [])
        if energies:
            min_e = min(energies)
            max_e = max(energies)
            barrier = max_e - min_e
            print(f"\n  Energy profile:")
            print(f"    Min energy:     {min_e:.4f} eV")
            print(f"    Max energy:     {max_e:.4f} eV")
            print(f"    Barrier:        {barrier:.4f} eV")
        if gaps:
            print(f"    Gap range:      {min(gaps):.4f} - {max(gaps):.4f} eV")

    if html_path.exists() and not args.no_open:
        print(f"\n[Launch] Opening reaction path animation...")
        try:
            webbrowser.open("file://" + str(html_path))
        except Exception:
            pass

    print("\n" + "=" * 64)


def main():
    parser = argparse.ArgumentParser(
        description="Hückel Molecular Orbital Calculator (Rust+Python+GPU)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "molecule", nargs="?", default=None,
        help="Path to molecule JSON file (atoms + bonds)"
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Output directory for results (default: results/<molecule_name>)"
    )
    parser.add_argument(
        "--no-gpu", action="store_true", default=False,
        help="Disable GPU acceleration (use CPU Jacobi method only)"
    )
    parser.add_argument(
        "--grid-res", type=int, default=50,
        help="Electron density grid resolution per axis (default: 50)"
    )
    parser.add_argument(
        "--grid-padding", type=float, default=3.0,
        help="Padding around molecule in Angstrom (default: 3.0)"
    )
    parser.add_argument(
        "--generate-graphene", metavar="NxN", default=None,
        help="Generate graphene sheet, e.g. 16x16 for 512 atoms"
    )
    parser.add_argument(
        "--no-open", action="store_true", default=False,
        help="Don't auto-open the HTML visualization"
    )
    parser.add_argument(
        "--rebuild", action="store_true", default=False,
        help="Force rebuild the Rust extension"
    )
    parser.add_argument(
        "--force-algorithm", choices=["auto", "simple", "ehmo"], default="auto",
        help="Force algorithm: 'auto' (detect automatically, default), 'simple' (pi-only), 'ehmo' (full valence + d-orbitals for metals)"
    )
    parser.add_argument(
        "--reaction", metavar="REACTION_JSON", default=None,
        help="Run reaction path simulation with a reaction JSON (reactant + product)"
    )
    parser.add_argument(
        "--num-steps", type=int, default=20,
        help="Number of IRC interpolation steps (default: 20, for --reaction mode)"
    )

    args = parser.parse_args()
    project_root = Path(__file__).parent.resolve()

    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

    print("=" * 64)
    print("  HUCKEL MOLECULAR ORBITAL CALCULATOR")
    print("  Rust Engine | wgpu GPU | Python Interface | Three.js Web")
    print("=" * 64)
    print()

    if args.rebuild or load_module(project_root) is None:
        print("[Setup] Need to build Rust extension.")
        if not check_rust_installed():
            print("[Error] Rust toolchain not detected.")
            print("        Install from: https://rustup.rs/")
            print("        Then re-run this script.")
            sys.exit(1)

        if not build_extension(project_root):
            print("[Error] Build failed.")
            sys.exit(1)

    module = load_module(project_root)
    if module is None:
        print("[Error] Failed to load extension after build.")
        sys.exit(1)

    if args.reaction:
        run_reaction_path_mode(args, module, project_root)
        return

    mol_path = None
    if args.generate_graphene:
        try:
            nx_str, ny_str = args.generate_graphene.lower().split("x")
            nx, ny = int(nx_str), int(ny_str)
            print(f"\n[Prep] Generating graphene sheet: {nx}x{ny} ({2*nx*ny} C atoms)")
            mol_path = generate_graphene(nx, ny, project_root)
        except Exception as e:
            print(f"[Error] Invalid graphene spec: {e}")
            print("        Use format like: --generate-graphene 16x16")
            sys.exit(1)
    elif args.molecule:
        mol_path = Path(args.molecule)
        if not mol_path.is_absolute():
            mol_path = (project_root / mol_path).resolve()
    else:
        print("Hint: Provide a molecule JSON file or use --generate-graphene")
        default_mol = project_root / "examples" / "molecules" / "benzene.json"
        if default_mol.exists():
            print(f"Using default: {default_mol}")
            mol_path = default_mol
        else:
            parser.print_help()
            sys.exit(0)

    if not mol_path.exists():
        print(f"[Error] Molecule file not found: {mol_path}")
        sys.exit(1)

    print(f"\n[Input] Molecule file: {mol_path}")

    try:
        with open(mol_path, "r", encoding="utf-8") as f:
            mol_data = json.load(f)
    except Exception as e:
        print(f"[Error] Cannot parse molecule JSON: {e}")
        sys.exit(1)

    user_force = args.force_algorithm
    algorithm, reason, heavy = detect_algorithm(mol_data, user_force)
    print(f"[Algorithm] {reason}")
    if heavy:
        print(f"[Algorithm] Heavy metal atoms detected: {', '.join(f'{e}(Z={z})' for e,z in heavy)}")
    print(f"[Algorithm] Final selection: {'Extended Hückel (EHMO)' if algorithm == 'ehmo' else 'Simple Hückel (pi-only)'}")

    mol_name = mol_path.stem
    output_dir = Path(args.output) if args.output else project_root / "results" / mol_name
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[Output] Results directory: {output_dir}")
    print(f"[Params] GPU acceleration: {not args.no_gpu}")
    print(f"[Params] Grid resolution: {args.grid_res}^3 (target)")
    print(f"[Params] Grid padding: {args.grid_padding} A")

    print("\n" + "-" * 64)
    t0 = time.time()

    try:
        mol_json_text = module.load_molecule_json(str(mol_path))

        engine_force = None
        if algorithm == "ehmo":
            engine_force = "ehmo"
        elif algorithm == "simple":
            engine_force = "simple"

        engine = module.HuckelEngine(use_gpu=not args.no_gpu, force_algorithm=engine_force)
        result = engine.calculate(
            mol_json_text,
            str(output_dir),
            grid_resolution=args.grid_res,
            grid_padding=args.grid_padding
        )
    except Exception as e:
        print(f"\n[Error] Calculation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    elapsed = time.time() - t0
    print("-" * 64)
    print(f"\n[Done] Total elapsed: {elapsed:.2f} seconds")

    html_path = output_dir / "visualization.html"
    svg_path = output_dir / "energy_levels.svg"
    bin_path = output_dir / "electron_density.bin"

    print("\nGenerated files:")
    for f in [svg_path, bin_path, html_path]:
        if f.exists():
            size_kb = f.stat().st_size / 1024
            print(f"  [OK] {f.relative_to(project_root)}  ({size_kb:.1f} KB)")

    if result:
        print("\nCalculation summary:")
        print(f"  Number of atoms:        {result.get('num_atoms', 'N/A')}")
        algo_used = result.get('algorithm', 'N/A')
        algo_display = {'SimpleHuckel': 'Simple Huckel (pi)', 'ExtendedHuckel': 'Extended Huckel (EHMO)'}.get(algo_used, algo_used)
        print(f"  Algorithm used:         {algo_display}")
        print(f"  Electrons:              {result.get('num_electrons', 'N/A')}")
        print(f"  HOMO-LUMO gap:          {result.get('homo_lumo_gap', 'N/A'):.4f} eV")
        evs = result.get('eigenvalues', [])
        if evs:
            print(f"  Eigenvalues (lowest 5): {[round(e, 4) for e in evs[:5]]}")
            if len(evs) > 5:
                print(f"  Eigenvalues (highest 3): {[round(e, 4) for e in evs[-3:]]}")
            import math
            has_imag = any(not math.isfinite(float(e)) for e in evs)
            if has_imag:
                print("[CRITICAL] Warning: non-finite values still in eigenvalues!")
            else:
                print("[OK] All eigenvalues are finite real numbers.")

    if html_path.exists() and not args.no_open:
        print(f"\n[Launch] Opening visualization in browser...")
        try:
            webbrowser.open("file://" + str(html_path))
        except Exception:
            pass

    print("\n" + "=" * 64)
    print("  Thank you for using Huckel Engine!")
    print("=" * 64)


if __name__ == "__main__":
    main()
