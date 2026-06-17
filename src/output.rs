use crate::types::{CalculationOutput, Molecule};
use std::fs::File;
use std::io::{Result, Write};
use byteorder::{LittleEndian, WriteBytesExt};
use base64::Engine;

pub fn write_density_binary(output: &CalculationOutput, path: &str) -> Result<()> {
    let mut file = File::create(path)?;

    let magic = b"DEN3";
    file.write_all(magic)?;

    file.write_u32::<LittleEndian>(1)?;

    file.write_u32::<LittleEndian>(output.grid_dims[0] as u32)?;
    file.write_u32::<LittleEndian>(output.grid_dims[1] as u32)?;
    file.write_u32::<LittleEndian>(output.grid_dims[2] as u32)?;

    file.write_f64::<LittleEndian>(output.grid_origin[0])?;
    file.write_f64::<LittleEndian>(output.grid_origin[1])?;
    file.write_f64::<LittleEndian>(output.grid_origin[2])?;

    file.write_f64::<LittleEndian>(output.grid_spacing)?;

    file.write_u32::<LittleEndian>(output.eigenvalues.len() as u32)?;
    for &e in &output.eigenvalues {
        file.write_f64::<LittleEndian>(e)?;
    }

    file.write_u32::<LittleEndian>(output.electron_density_grid.len() as u32)?;
    for &d in &output.electron_density_grid {
        file.write_f32::<LittleEndian>(d)?;
    }

    Ok(())
}

pub fn density_to_base64(output: &CalculationOutput) -> String {
    let mut buf = Vec::new();

    let _ = buf.write_u32::<LittleEndian>(output.grid_dims[0] as u32);
    let _ = buf.write_u32::<LittleEndian>(output.grid_dims[1] as u32);
    let _ = buf.write_u32::<LittleEndian>(output.grid_dims[2] as u32);

    let _ = buf.write_f64::<LittleEndian>(output.grid_origin[0]);
    let _ = buf.write_f64::<LittleEndian>(output.grid_origin[1]);
    let _ = buf.write_f64::<LittleEndian>(output.grid_origin[2]);

    let _ = buf.write_f64::<LittleEndian>(output.grid_spacing);

    for &d in &output.electron_density_grid {
        let _ = buf.write_f32::<LittleEndian>(d);
    }

    base64::engine::general_purpose::STANDARD.encode(&buf)
}

pub fn molecule_to_json(mol: &Molecule, output: &CalculationOutput) -> String {
    let atoms: Vec<String> = mol.atoms.iter().map(|a| {
        format!(
            r#"{{"element":"{}","x":{:.6},"y":{:.6},"z":{:.6},"radius":{:.3}}}"#,
            a.element, a.x, a.y, a.z,
            Molecule::atom_radius(&a.element)
        )
    }).collect();

    let eigenvalues_str: Vec<String> = output.eigenvalues.iter()
        .map(|e| format!("{:.6}", e)).collect();

    format!(
        r#"{{"atoms":[{}],"grid_dims":[{},{},{}],"origin":[{:.4},{:.4},{:.4}],"spacing":{:.4},"eigenvalues":[{}],"num_electrons":{},"homo_lumo_gap":{:.6}}}"#,
        atoms.join(","),
        output.grid_dims[0], output.grid_dims[1], output.grid_dims[2],
        output.grid_origin[0], output.grid_origin[1], output.grid_origin[2],
        output.grid_spacing,
        eigenvalues_str.join(","),
        output.num_electrons,
        output.homo_lumo_gap
    )
}
