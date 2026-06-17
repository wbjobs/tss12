pub fn get_html_template(
    molecule_json: &str,
    density_b64: &str,
    svg_content: &str,
) -> String {
    let svg_b64 = base64::engine::general_purpose::STANDARD.encode(svg_content.as_bytes());
    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>Hückel Quantum Mechanics Visualization</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background: linear-gradient(135deg, #0f0c29, #302b63, #24243e);
    color: #e0e0e0;
    min-height: 100vh;
    overflow-x: hidden;
  }}
  .header {{
    padding: 24px 32px;
    background: rgba(0,0,0,0.3);
    border-bottom: 1px solid rgba(255,255,255,0.1);
  }}
  .header h1 {{
    font-size: 24px;
    background: linear-gradient(90deg, #60a5fa, #a78bfa, #f472b6);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    margin-bottom: 8px;
  }}
  .header p {{ font-size: 13px; color: #9ca3af; }}
  .container {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 20px;
    padding: 20px;
  }}
  @media (max-width: 1200px) {{
    .container {{ grid-template-columns: 1fr; }}
  }}
  .panel {{
    background: rgba(255,255,255,0.04);
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 12px;
    padding: 20px;
    backdrop-filter: blur(10px);
  }}
  .panel-title {{
    font-size: 16px;
    font-weight: 600;
    color: #c4b5fd;
    margin-bottom: 16px;
    padding-bottom: 10px;
    border-bottom: 1px solid rgba(196,181,253,0.2);
    display: flex;
    align-items: center;
    gap: 8px;
  }}
  .panel-title::before {{
    content: '';
    width: 4px;
    height: 16px;
    background: linear-gradient(180deg, #60a5fa, #a78bfa);
    border-radius: 2px;
  }}
  #viewer3d {{
    width: 100%;
    height: 520px;
    border-radius: 8px;
    background: radial-gradient(ellipse at center, #1a1a2e 0%, #0a0a14 100%);
  }}
  #energyLevelDiagram {{
    width: 100%;
    height: 520px;
    overflow: auto;
    border-radius: 8px;
    background: white;
  }}
  #energyLevelDiagram svg {{ width: 100%; height: auto; min-height: 100%; }}
  .controls {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin-top: 16px;
  }}
  .control-group {{
    display: flex;
    flex-direction: column;
    gap: 6px;
    flex: 1;
    min-width: 140px;
  }}
  .control-group label {{
    font-size: 12px;
    color: #9ca3af;
  }}
  .control-group input[type=range] {{
    width: 100%;
    accent-color: #a78bfa;
  }}
  .control-group .value-display {{
    font-size: 11px;
    color: #c4b5fd;
    text-align: right;
    font-family: monospace;
  }}
  .checkbox-group {{
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 13px;
    color: #d1d5db;
    cursor: pointer;
    padding: 6px 0;
  }}
  .checkbox-group input {{ accent-color: #a78bfa; width: 16px; height: 16px; }}
  .info-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 12px;
    margin-top: 16px;
  }}
  .info-card {{
    background: rgba(167,139,250,0.08);
    border: 1px solid rgba(167,139,250,0.15);
    border-radius: 8px;
    padding: 12px;
  }}
  .info-label {{ font-size: 11px; color: #9ca3af; text-transform: uppercase; letter-spacing: 0.5px; }}
  .info-value {{ font-size: 18px; font-weight: 700; color: #e0e7ff; margin-top: 4px; font-family: 'SF Mono', monospace; }}
  .loading {{
    display: flex;
    align-items: center;
    justify-content: center;
    height: 100%;
    color: #a78bfa;
    font-size: 14px;
  }}
  .spinner {{
    border: 2px solid rgba(167,139,250,0.2);
    border-top-color: #a78bfa;
    border-radius: 50%;
    width: 24px;
    height: 24px;
    animation: spin 0.8s linear infinite;
    margin-right: 10px;
  }}
  @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
</style>
</head>
<body>
<div class="header">
  <h1>⚛️ Hückel Quantum Mechanics Calculation Results</h1>
  <p>Electron density isosurface visualization & molecular orbital energy level diagram</p>
</div>

<div class="container">
  <div class="panel">
    <div class="panel-title">3D Electron Density Isosurface</div>
    <div id="viewer3d">
      <div class="loading" id="loading3d">
        <div class="spinner"></div>Initializing WebGL renderer...
      </div>
    </div>
    <div class="controls">
      <div class="control-group">
        <label>Isosurface Level</label>
        <input type="range" id="isoLevel" min="0.005" max="0.3" step="0.005" value="0.05"/>
        <div class="value-display" id="isoLevelValue">0.050</div>
      </div>
      <div class="control-group">
        <label>Opacity</label>
        <input type="range" id="opacity" min="0.1" max="1.0" step="0.05" value="0.75"/>
        <div class="value-display" id="opacityValue">0.75</div>
      </div>
      <div class="control-group">
        <label>View Mode</label>
        <div style="display:flex;gap:6px;flex-wrap:wrap;">
          <label class="checkbox-group"><input type="checkbox" id="showDensity" checked/>Density</label>
          <label class="checkbox-group"><input type="checkbox" id="showAtoms" checked/>Atoms</label>
          <label class="checkbox-group"><input type="checkbox" id="showBonds" checked/>Bonds</label>
        </div>
      </div>
    </div>
    <div class="info-grid" id="infoGrid"></div>
  </div>

  <div class="panel">
    <div class="panel-title">Molecular Orbital Energy Levels</div>
    <div id="energyLevelDiagram"></div>
  </div>
</div>

<script type="importmap">
{{
  "imports": {{
    "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"
  }}
}}
</script>

<script type="module">
import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';
import {{ MarchingCubes }} from 'three/addons/objects/MarchingCubes.js';

const MOLECULE_DATA = {molecule_json};
const DENSITY_B64 = "{density_b64}";
const SVG_B64 = "{svg_b64}";

const ELEMENT_COLORS = {{
  'H': 0xffffff, 'C': 0x303030, 'N': 0x3050f8, 'O': 0xff0d0d,
  'F': 0x90e050, 'P': 0xff8000, 'S': 0xffff30, 'Cl': 0x1ff01f,
  'B': 0xffb5b5, 'Br': 0xa62929, 'I': 0x940094
}};

function b64ToBytes(b64) {{
  const binary = atob(b64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return bytes;
}}

function initEnergyDiagram() {{
  const svg = atob(SVG_B64);
  document.getElementById('energyLevelDiagram').innerHTML = svg;
}}

function parseDensityData() {{
  const bytes = b64ToBytes(DENSITY_B64);
  const dv = new DataView(bytes.buffer);
  let offset = 0;
  const dims = [dv.getUint32(offset, true), dv.getUint32(offset+4, true), dv.getUint32(offset+8, true)];
  offset += 12;
  const origin = [dv.getFloat64(offset, true), dv.getFloat64(offset+8, true), dv.getFloat64(offset+16, true)];
  offset += 24;
  const spacing = dv.getFloat64(offset, true);
  offset += 8;
  const nVoxels = dims[0]*dims[1]*dims[2];
  const density = new Float32Array(nVoxels);
  for (let i = 0; i < nVoxels; i++) {{
    density[i] = dv.getFloat32(offset, true);
    offset += 4;
  }}
  return {{ dims, origin, spacing, density }};
}}

function initInfoPanel() {{
  const grid = document.getElementById('infoGrid');
  const eigenvalues = MOLECULE_DATA.eigenvalues;
  const n = MOLECULE_DATA.num_electrons;
  const homoIdx = Math.floor((n-1)/2);
  const lumoIdx = homoIdx + 1;
  const gap = MOLECULE_DATA.homo_lumo_gap;
  const cards = [
    {{ label: 'Atoms', value: MOLECULE_DATA.atoms.length }},
    {{ label: 'π-Electrons', value: n }},
    {{ label: 'Orbitals', value: eigenvalues.length }},
    {{ label: 'HOMO (eV)', value: eigenvalues[homoIdx]?.toFixed(4) || 'N/A' }},
    {{ label: 'LUMO (eV)', value: eigenvalues[lumoIdx]?.toFixed(4) || 'N/A' }},
    {{ label: 'HOMO-LUMO Gap', value: gap.toFixed(4) + ' eV' }},
    {{ label: 'Grid Dim X', value: MOLECULE_DATA.grid_dims[0] }},
    {{ label: 'Grid Dim Y', value: MOLECULE_DATA.grid_dims[1] }},
    {{ label: 'Grid Dim Z', value: MOLECULE_DATA.grid_dims[2] }},
  ];
  grid.innerHTML = cards.map(c => `
    <div class="info-card">
      <div class="info-label">${{c.label}}</div>
      <div class="info-value">${{c.value}}</div>
    </div>
  `).join('');
}}

function initThreeJS() {{
  const container = document.getElementById('viewer3d');
  const loading = document.getElementById('loading3d');

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0a0a14);
  scene.fog = new THREE.FogExp2(0x0a0a14, 0.015);

  const camera = new THREE.PerspectiveCamera(60, container.clientWidth/container.clientHeight, 0.1, 1000);
  camera.position.set(8, 6, 12);

  const renderer = new THREE.WebGLRenderer({{ antialias: true, alpha: true }});
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  container.innerHTML = '';
  container.appendChild(renderer.domElement);

  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;

  scene.add(new THREE.AmbientLight(0x404050, 0.6));
  const keyLight = new THREE.DirectionalLight(0xffffff, 1.2);
  keyLight.position.set(10, 15, 8);
  keyLight.castShadow = true;
  keyLight.shadow.mapSize.width = 2048;
  keyLight.shadow.mapSize.height = 2048;
  scene.add(keyLight);
  const fillLight = new THREE.DirectionalLight(0x8080ff, 0.5);
  fillLight.position.set(-8, 5, -6);
  scene.add(fillLight);
  const rimLight = new THREE.PointLight(0xff80ff, 0.6, 50);
  rimLight.position.set(-5, -3, -8);
  scene.add(rimLight);

  const densityData = parseDensityData();
  const {{ dims, origin, spacing, density }} = densityData;

  let maxDensity = 0;
  for (let i = 0; i < density.length; i++) {{
    if (density[i] > maxDensity) maxDensity = density[i];
  }}
  console.log('[3D] Max density:', maxDensity.toFixed(4));

  function getDensity(x, y, z) {{
    const ix = Math.floor((x - origin[0]) / spacing);
    const iy = Math.floor((y - origin[1]) / spacing);
    const iz = Math.floor((z - origin[2]) / spacing);
    if (ix < 0 || ix >= dims[0] || iy < 0 || iy >= dims[1] || iz < 0 || iz >= dims[2]) return 0;
    return density[iz * dims[1] * dims[0] + iy * dims[0] + ix];
  }}

  const resolution = 80;
  const effect = new MarchingCubes(resolution, null, true, true, 500000);
  effect.enableUvs = false;
  effect.enableColors = true;
  const isoMaterial = new THREE.MeshPhongMaterial({{
    color: 0x60a5fa,
    transparent: true,
    opacity: 0.75,
    shininess: 80,
    specular: 0x404080,
    side: THREE.DoubleSide,
    vertexColors: true,
  }});
  effect.material = isoMaterial;

  const size = {{
    x: (dims[0] - 1) * spacing,
    y: (dims[1] - 1) * spacing,
    z: (dims[2] - 1) * spacing,
  }};
  effect.scale.set(size.x, size.y, size.z);
  effect.position.set(origin[0] + size.x/2, origin[1] + size.y/2, origin[2] + size.z/2);

  function updateIsosurface(level) {{
    effect.reset();
    const field = new Float32Array(resolution ** 3);
    for (let iz = 0; iz < resolution; iz++) {{
      const z = origin[2] + (iz / (resolution - 1)) * size.z;
      for (let iy = 0; iy < resolution; iy++) {{
        const y = origin[1] + (iy / (resolution - 1)) * size.y;
        for (let ix = 0; ix < resolution; ix++) {{
          const x = origin[0] + (ix / (resolution - 1)) * size.x;
          const d = getDensity(x, y, z);
          field[iz*resolution*resolution + iy*resolution + ix] = d - level;
        }}
      }}
    }}
    effect.generate(field, level);
    const geometry = effect.geometry;
    if (geometry.attributes.position && geometry.attributes.position.count > 0) {{
      const colors = new Float32Array(geometry.attributes.position.count * 3);
      for (let i = 0; i < geometry.attributes.position.count; i++) {{
        colors[i*3]   = 0.376 + 0.2 * Math.random();
        colors[i*3+1] = 0.647 + 0.15 * Math.random();
        colors[i*3+2] = 0.980 + 0.02 * Math.random();
      }}
      geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    }}
  }}

  const moleculeGroup = new THREE.Group();
  const atoms = MOLECULE_DATA.atoms;

  atoms.forEach(atom => {{
    const color = ELEMENT_COLORS[atom.element] || 0xcccccc;
    const radius = Math.max(0.2, atom.radius * 0.5);
    const geom = new THREE.SphereGeometry(radius, 32, 32);
    const mat = new THREE.MeshPhongMaterial({{
      color, shininess: 120, specular: 0x888888
    }});
    const sphere = new THREE.Mesh(geom, mat);
    sphere.position.set(atom.x, atom.y, atom.z);
    sphere.castShadow = true;
    sphere.receiveShadow = true;
    moleculeGroup.add(sphere);
  }});

  for (let i = 0; i < atoms.length; i++) {{
    for (let j = i+1; j < atoms.length; j++) {{
      const a = atoms[i], b = atoms[j];
      const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
      const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
      const cov = a.radius + b.radius;
      if (dist < cov * 1.2) {{
        const dir = new THREE.Vector3(dx, dy, dz);
        const bondGeom = new THREE.CylinderGeometry(0.06, 0.06, dist, 12);
        const bondMat = new THREE.MeshPhongMaterial({{ color: 0x888899, shininess: 100 }});
        const bond = new THREE.Mesh(bondGeom, bondMat);
        bond.position.set((a.x+b.x)/2, (a.y+b.y)/2, (a.z+b.z)/2);
        const up = new THREE.Vector3(0, 1, 0);
        const q = new THREE.Quaternion().setFromUnitVectors(up, dir.clone().normalize());
        bond.quaternion.copy(q);
        bond.castShadow = true;
        moleculeGroup.add(bond);
      }}
    }}
  }}

  scene.add(moleculeGroup);
  scene.add(effect);

  const box = new THREE.Box3().setFromObject(moleculeGroup);
  const center = box.getCenter(new THREE.Vector3());
  const radius = box.getSize(new THREE.Vector3()).length() * 0.7;
  controls.target.copy(center);
  camera.position.copy(center).add(new THREE.Vector3(radius, radius*0.7, radius*1.2));
  controls.update();

  let currentLevel = 0.05;
  updateIsosurface(currentLevel);

  const isoSlider = document.getElementById('isoLevel');
  const isoValue = document.getElementById('isoLevelValue');
  const opacitySlider = document.getElementById('opacity');
  const opacityValue = document.getElementById('opacityValue');
  const showDensity = document.getElementById('showDensity');
  const showAtoms = document.getElementById('showAtoms');
  const showBonds = document.getElementById('showBonds');

  let updateTimer = null;
  isoSlider.addEventListener('input', (e) => {{
    currentLevel = parseFloat(e.target.value);
    isoValue.textContent = currentLevel.toFixed(3);
    if (updateTimer) clearTimeout(updateTimer);
    updateTimer = setTimeout(() => updateIsosurface(currentLevel), 100);
  }});
  opacitySlider.addEventListener('input', (e) => {{
    const val = parseFloat(e.target.value);
    opacityValue.textContent = val.toFixed(2);
    isoMaterial.opacity = val;
  }});
  showDensity.addEventListener('change', () => {{ effect.visible = showDensity.checked; }});
  showAtoms.addEventListener('change', () => {{
    moleculeGroup.children.forEach(c => {{
      if (c.geometry instanceof THREE.SphereGeometry) c.visible = showAtoms.checked;
    }});
  }});
  showBonds.addEventListener('change', () => {{
    moleculeGroup.children.forEach(c => {{
      if (c.geometry instanceof THREE.CylinderGeometry) c.visible = showBonds.checked;
    }});
  }});

  window.addEventListener('resize', () => {{
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  }});

  function animate() {{
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
  }}
  animate();
}}

initEnergyDiagram();
initInfoPanel();
initThreeJS();
</script>
</body>
</html>"#
    )
}
