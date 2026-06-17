pub fn generate_energy_level_svg(
    eigenvalues: &Vec<f64>,
    num_electrons: usize,
    width: usize,
    height: usize,
) -> String {
    let n = eigenvalues.len();
    let homo_index = if num_electrons > 0 {
        (num_electrons - 1) / 2
    } else {
        0
    };
    let lumo_index = homo_index + 1;

    let min_e = eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_e = eigenvalues.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let e_range = (max_e - min_e).max(0.1);
    let padding = 60.0;
    let plot_width = width as f64 - 2.0 * padding;
    let plot_height = height as f64 - 2.0 * padding;
    let line_width = plot_width * 0.3;

    let mut svg = String::new();
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">"#,
        width, height, width, height
    ));
    svg.push_str(r#"<defs><style>"#);
    svg.push_str(".title { font: bold 18px sans-serif; fill: #1a1a2e; }");
    svg.push_str(".axis-label { font: 12px sans-serif; fill: #333; }");
    svg.push_str(".homo-line { stroke: #3b82f6; stroke-width: 3; }");
    svg.push_str(".lumo-line { stroke: #ef4444; stroke-width: 3; stroke-dasharray: 8,4; }");
    svg.push_str(".level-line { stroke: #64748b; stroke-width: 2.5; }");
    svg.push_str(".grid-line { stroke: #e2e8f0; stroke-width: 1; }");
    svg.push_str(".gap-label { font: bold 14px sans-serif; fill: #7c3aed; }");
    svg.push_str("</style></defs>");

    svg.push_str(&format!(
        r#"<rect width="100%" height="100%" fill="#f8fafc"/>"#
    ));

    svg.push_str(&format!(
        r#"<text x="{}" y="30" text-anchor="middle" class="title">Hückel Molecular Orbital Energy Levels</text>"#,
        width as f64 / 2.0
    ));

    svg.push_str(&format!(
        r#"<text x="20" y="{}" transform="rotate(-90 20 {})" text-anchor="middle" class="axis-label">Energy (eV)</text>"#,
        height as f64 / 2.0,
        height as f64 / 2.0
    ));

    let num_ticks = 8;
    for t in 0..=num_ticks {
        let frac = t as f64 / num_ticks as f64;
        let y = padding + frac * plot_height;
        let e_val = min_e + (1.0 - frac) * e_range;
        svg.push_str(&format!(
            r#"<line x1="{}" y1="{}" x2="{}" y2="{}" class="grid-line"/>"#,
            padding - 10.0, y, padding + plot_width + 10.0, y
        ));
        svg.push_str(&format!(
            r#"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" class="axis-label">{:.2}</text>"#,
            padding - 15.0, y, e_val
        ));
    }

    let center_x = padding + plot_width / 2.0;
    for (i, &e) in eigenvalues.iter().enumerate() {
        let frac = (e - min_e) / e_range;
        let y = padding + (1.0 - frac) * plot_height;

        let is_homo = i == homo_index && i < n;
        let is_lumo = i == lumo_index && i < n;

        let (x1, x2, class) = if is_homo {
            (
                center_x - line_width / 2.0,
                center_x + line_width / 2.0,
                "homo-line",
            )
        } else if is_lumo {
            (
                center_x - line_width / 2.0,
                center_x + line_width / 2.0,
                "lumo-line",
            )
        } else {
            (
                center_x - line_width / 2.0,
                center_x + line_width / 2.0,
                "level-line",
            )
        };

        svg.push_str(&format!(
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" class="{}"/>"#,
            x1, y, x2, y, class
        ));

        if is_homo || is_lumo {
            let label = if is_homo { "HOMO" } else { "LUMO" };
            let label_color = if is_homo { "#1e40af" } else { "#991b1b" };
            svg.push_str(&format!(
                r#"<text x="{:.2}" y="{:.2}" text-anchor="start" dominant-baseline="middle" 
                      font-size="12" font-weight="bold" fill="{}" dx="12">{}</text>"#,
                x2, y, label_color, label
            ));
            svg.push_str(&format!(
                r#"<text x="{:.2}" y="{:.2}" text-anchor="end" dominant-baseline="middle" 
                      font-size="11" fill="#555" dx="-12">{:.4} eV</text>"#,
                x1, y, e
            ));
        }

        if num_electrons > 0 && i < homo_index {
            let arrow_y1 = y - 8.0;
            let arrow_y2 = y + 8.0;
            svg.push_str(&format!(
                r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" 
                      stroke="#10b981" stroke-width="1.5" marker-end="url(#uparrow)"/>"#,
                center_x - 5.0, arrow_y2, center_x - 5.0, arrow_y1
            ));
            svg.push_str(&format!(
                r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" 
                      stroke="#10b981" stroke-width="1.5"/>"#,
                center_x + 5.0, arrow_y1, center_x + 5.0, arrow_y2
            ));
        } else if num_electrons > 0 && i == homo_index {
            let arrow_y1 = y - 8.0;
            let arrow_y2 = y + 8.0;
            svg.push_str(&format!(
                r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" 
                      stroke="#10b981" stroke-width="1.5"/>"#,
                center_x, arrow_y2, center_x, arrow_y1
            ));
            if num_electrons % 2 == 0 {
                svg.push_str(&format!(
                    r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" 
                          stroke="#10b981" stroke-width="1.5"/>"#,
                    center_x + 10.0, arrow_y1, center_x + 10.0, arrow_y2
                ));
            }
        }
    }

    if lumo_index < n && homo_index < n {
        let gap = eigenvalues[lumo_index] - eigenvalues[homo_index];
        let homo_frac = (eigenvalues[homo_index] - min_e) / e_range;
        let lumo_frac = (eigenvalues[lumo_index] - min_e) / e_range;
        let y_homo = padding + (1.0 - homo_frac) * plot_height;
        let y_lumo = padding + (1.0 - lumo_frac) * plot_height;

        svg.push_str(&format!(
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" 
                  stroke="#7c3aed" stroke-width="1" stroke-dasharray="4,4"/>"#,
            center_x + line_width / 2.0 + 50.0, y_homo,
            center_x + line_width / 2.0 + 50.0, y_lumo
        ));
        svg.push_str(&format!(
            r#"<text x="{:.2}" y="{:.2}" text-anchor="middle" dominant-baseline="middle" 
                  class="gap-label" transform="rotate(90 {:.2} {:.2})">
                  ΔE (HOMO-LUMO) = {:.4} eV</text>"#,
            center_x + line_width / 2.0 + 70.0,
            (y_homo + y_lumo) / 2.0,
            center_x + line_width / 2.0 + 70.0,
            (y_homo + y_lumo) / 2.0,
            gap
        ));
    }

    svg.push_str(&format!(
        r#"<text x="{}" y="{}" text-anchor="middle" class="axis-label">
              N = {} atoms, N_e = {} π-electrons</text>"#,
        width as f64 / 2.0,
        height as f64 - 20.0,
        n,
        num_electrons
    ));

    svg.push_str("</svg>");
    svg
}
