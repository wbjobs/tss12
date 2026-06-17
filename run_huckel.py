#!/usr/bin/env python3
"""
Hückel Quantum Mechanics Calculation Runner

通过PyO3调用Rust库执行半经验量子力学计算。

用法:
    python run_huckel.py <molecule.json> [--output DIR] [--no-gpu]
                         [--grid-res N] [--grid-padding Å]
                         [--generate-graphene NxN]

示例:
    python run_huckel.py examples/molecules/benzene.json
    python run_huckel.py examples/molecules/naphthalene.json --output results/naphthalene
    python run_huckel.py --generate-graphene 16x16 --grid-res 64
    python run_huckel.py examples/molecules/benzene.json --no-gpu
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

    py_ver = f"cp{sys.version_info.major}{sys.version_info.minor}"
    abi_tag = f"{py_ver}-{py_ver}"
    plat_tag = {
        "win32": "win_amd64",
        "darwin": "macosx_11_0_arm64" if "arm" in os.uname().machine.lower() else "macosx_10_9_x86_64",
        "linux": "linux_x86_64",
    }.get(sys.platform, "any")

    src_lib = project_root / "target" / "release" / lib_name
    if not src_lib.exists():
        print(f"[Error] Built library not found at {src_lib}")
        return False

    import shutil
    import platform

    ext_suffix = ""
    try:
        import sysconfig
        ext_suffix = sysconfig.get_config_var("EXT_SUFFIX") or f".{lib_ext}"
    except:
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
        help="Padding around molecule in Å (default: 3.0)"
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

    args = parser.parse_args()
    project_root = Path(__file__).parent.resolve()

    print("=" * 64)
    print("  HÜCKEL MOLECULAR ORBITAL CALCULATOR")
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

    mol_name = mol_path.stem
    output_dir = Path(args.output) if args.output else project_root / "results" / mol_name
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[Input] Molecule file: {mol_path}")
    print(f"[Output] Results directory: {output_dir}")
    print(f"[Params] GPU acceleration: {not args.no_gpu}")
    print(f"[Params] Grid resolution: {args.grid_res}^3 (target)")
    print(f"[Params] Grid padding: {args.grid_padding} Å")

    print("\n" + "-" * 64)
    t0 = time.time()

    try:
        mol_json_text = module.load_molecule_json(str(mol_path))

        engine = module.HuckelEngine(use_gpu=not args.no_gpu)
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
            print(f"  ✓ {f.relative_to(project_root)}  ({size_kb:.1f} KB)")

    if result:
        print("\nCalculation summary:")
        print(f"  Number of atoms:        {result.get('num_atoms', 'N/A')}")
        print(f"  π-electrons:            {result.get('num_electrons', 'N/A')}")
        print(f"  HOMO-LUMO gap:          {result.get('homo_lumo_gap', 'N/A'):.4f} eV")
        evs = result.get('eigenvalues', [])
        if evs:
            print(f"  Eigenvalues (lowest 5): {[round(e, 4) for e in evs[:5]]}")
            if len(evs) > 5:
                print(f"  Eigenvalues (highest 3): {[round(e, 4) for e in evs[-3:]]}")

    if html_path.exists() and not args.no_open:
        print(f"\n[Launch] Opening visualization in browser...")
        try:
            webbrowser.open("file://" + str(html_path))
        except Exception:
            pass

    print("\n" + "=" * 64)
    print("  Thank you for using Hückel Engine!")
    print("=" * 64)


if __name__ == "__main__":
    main()
