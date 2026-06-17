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
  @media (max-width: 1200px) {{ .container {{ grid-template-columns: 1fr; }} }}
  .panel {{
    background: rgba(255,255,255,0.04);
    border: 1px solid rgba(255,255,255,0.08);
    border-radius: 12px;
    padding: 20px;
    backdrop-filter: blur(10px);
  }}
  .panel-title {{
    font-size: 16px; font-weight: 600; color: #c4b5fd;
    margin-bottom: 16px; padding-bottom: 10px;
    border-bottom: 1px solid rgba(196,181,253,0.2);
    display: flex; align-items: center; gap: 8px;
  }}
  .panel-title::before {{
    content: ''; width: 4px; height: 16px;
    background: linear-gradient(180deg, #60a5fa, #a78bfa); border-radius: 2px;
  }}
  #viewer3d {{ width:100%; height:520px; border-radius:8px; background:radial-gradient(ellipse at center,#1a1a2e 0%,#0a0a14 100%); }}
  #energyLevelDiagram {{ width:100%; height:520px; overflow:auto; border-radius:8px; background:white; }}
  #energyLevelDiagram svg {{ width:100%; height:auto; min-height:100%; }}
  .controls {{ display:flex; flex-wrap:wrap; gap:10px; margin-top:16px; }}
  .control-group {{ display:flex; flex-direction:column; gap:6px; flex:1; min-width:140px; }}
  .control-group label {{ font-size:12px; color:#9ca3af; }}
  .control-group input[type=range] {{ width:100%; accent-color:#a78bfa; }}
  .control-group .value-display {{ font-size:11px; color:#c4b5fd; text-align:right; font-family:monospace; }}
  .checkbox-group {{ display:flex; align-items:center; gap:8px; font-size:13px; color:#d1d5db; cursor:pointer; padding:6px 0; }}
  .checkbox-group input {{ accent-color:#a78bfa; width:16px; height:16px; }}
  .info-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); gap:12px; margin-top:16px; }}
  .info-card {{ background:rgba(167,139,250,0.08); border:1px solid rgba(167,139,250,0.15); border-radius:8px; padding:12px; }}
  .info-label {{ font-size:11px; color:#9ca3af; text-transform:uppercase; letter-spacing:0.5px; }}
  .info-value {{ font-size:18px; font-weight:700; color:#e0e7ff; margin-top:4px; font-family:'SF Mono',monospace; }}
  .loading {{ display:flex; align-items:center; justify-content:center; height:100%; color:#a78bfa; font-size:14px; }}
  .spinner {{ border:2px solid rgba(167,139,250,0.2); border-top-color:#a78bfa; border-radius:50%; width:24px; height:24px; animation:spin 0.8s linear infinite; margin-right:10px; }}
  @keyframes spin {{ to {{ transform:rotate(360deg); }} }}
</style>
</head>
<body>
<div class="header">
  <h1>Huckel Quantum Mechanics Calculation Results</h1>
  <p>Electron density isosurface visualization & molecular orbital energy level diagram</p>
</div>
<div class="container">
  <div class="panel">
    <div class="panel-title">3D Electron Density Isosurface</div>
    <div id="viewer3d"><div class="loading" id="loading3d"><div class="spinner"></div>Initializing WebGL...</div></div>
    <div class="controls">
      <div class="control-group"><label>Isosurface Level</label><input type="range" id="isoLevel" min="0.005" max="0.3" step="0.005" value="0.05"/><div class="value-display" id="isoLevelValue">0.050</div></div>
      <div class="control-group"><label>Opacity</label><input type="range" id="opacity" min="0.1" max="1.0" step="0.05" value="0.75"/><div class="value-display" id="opacityValue">0.75</div></div>
      <div class="control-group"><label>View Mode</label><div style="display:flex;gap:6px;flex-wrap:wrap;">
        <label class="checkbox-group"><input type="checkbox" id="showDensity" checked/>Density</label>
        <label class="checkbox-group"><input type="checkbox" id="showAtoms" checked/>Atoms</label>
        <label class="checkbox-group"><input type="checkbox" id="showBonds" checked/>Bonds</label>
      </div></div>
    </div>
    <div class="info-grid" id="infoGrid"></div>
  </div>
  <div class="panel">
    <div class="panel-title">Molecular Orbital Energy Levels</div>
    <div id="energyLevelDiagram"></div>
  </div>
</div>
<script type="importmap">{{ "imports": {{ "three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/" }} }}</script>
<script type="module">
import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';
import {{ MarchingCubes }} from 'three/addons/objects/MarchingCubes.js';
const MOLECULE_DATA = {molecule_json};
const DENSITY_B64 = "{density_b64}";
const SVG_B64 = "{svg_b64}";
const ELEMENT_COLORS = {{ 'H':0xffffff,'C':0x303030,'N':0x3050f8,'O':0xff0d0d,'F':0x90e050,'P':0xff8000,'S':0xffff30,'Cl':0x1ff01f,'B':0xffb5b5,'Br':0xa62929,'I':0x940094 }};
function b64ToBytes(b64){{ const binary=atob(b64); const bytes=new Uint8Array(binary.length); for(let i=0;i<binary.length;i++) bytes[i]=binary.charCodeAt(i); return bytes; }}
function initEnergyDiagram(){{ document.getElementById('energyLevelDiagram').innerHTML=atob(SVG_B64); }}
function parseDensityData(){{ const bytes=b64ToBytes(DENSITY_B64); const dv=new DataView(bytes.buffer); let offset=0; const dims=[dv.getUint32(offset,true),dv.getUint32(offset+4,true),dv.getUint32(offset+8,true)]; offset+=12; const origin=[dv.getFloat64(offset,true),dv.getFloat64(offset+8,true),dv.getFloat64(offset+16,true)]; offset+=24; const spacing=dv.getFloat64(offset,true); offset+=8; const nVoxels=dims[0]*dims[1]*dims[2]; const density=new Float32Array(nVoxels); for(let i=0;i<nVoxels;i++){{ density[i]=dv.getFloat32(offset,true); offset+=4; }} return {{dims,origin,spacing,density}}; }}
function initInfoPanel(){{
  const grid=document.getElementById('infoGrid');
  const eigenvalues=MOLECULE_DATA.eigenvalues; const n=MOLECULE_DATA.num_electrons;
  const homoIdx=Math.floor((n-1)/2); const lumoIdx=homoIdx+1;
  const cards=[
    {{label:'Atoms',value:MOLECULE_DATA.atoms.length}},
    {{label:'Electrons',value:n}},
    {{label:'Orbitals',value:eigenvalues.length}},
    {{label:'HOMO (eV)',value:eigenvalues[homoIdx]?.toFixed(4)||'N/A'}},
    {{label:'LUMO (eV)',value:eigenvalues[lumoIdx]?.toFixed(4)||'N/A'}},
    {{label:'HOMO-LUMO Gap',value:MOLECULE_DATA.homo_lumo_gap.toFixed(4)+' eV'}},
  ];
  grid.innerHTML=cards.map(c=>`<div class="info-card"><div class="info-label">${{c.label}}</div><div class="info-value">${{c.value}}</div></div>`).join('');
}}
function initThreeJS(){{
  const container=document.getElementById('viewer3d');
  const scene=new THREE.Scene(); scene.background=new THREE.Color(0x0a0a14); scene.fog=new THREE.FogExp2(0x0a0a14,0.015);
  const camera=new THREE.PerspectiveCamera(60,container.clientWidth/container.clientHeight,0.1,1000); camera.position.set(8,6,12);
  const renderer=new THREE.WebGLRenderer({{antialias:true,alpha:true}}); renderer.setSize(container.clientWidth,container.clientHeight); renderer.setPixelRatio(Math.min(window.devicePixelRatio,2)); renderer.shadowMap.enabled=true; renderer.shadowMap.type=THREE.PCFSoftShadowMap;
  container.innerHTML=''; container.appendChild(renderer.domElement);
  const controls=new OrbitControls(camera,renderer.domElement); controls.enableDamping=true; controls.dampingFactor=0.08;
  scene.add(new THREE.AmbientLight(0x404050,0.6));
  const keyLight=new THREE.DirectionalLight(0xffffff,1.2); keyLight.position.set(10,15,8); keyLight.castShadow=true; keyLight.shadow.mapSize.width=2048; keyLight.shadow.mapSize.height=2048; scene.add(keyLight);
  scene.add(new THREE.DirectionalLight(0x8080ff,0.5).position.set(-8,5,-6)||new THREE.DirectionalLight(0x8080ff,0.5));
  const rimLight=new THREE.PointLight(0xff80ff,0.6,50); rimLight.position.set(-5,-3,-8); scene.add(rimLight);
  const densityData=parseDensityData(); const {{dims,origin,spacing,density}}=densityData;
  let maxDensity=0; for(let i=0;i<density.length;i++) if(density[i]>maxDensity) maxDensity=density[i];
  function getDensity(x,y,z){{ const ix=Math.floor((x-origin[0])/spacing); const iy=Math.floor((y-origin[1])/spacing); const iz=Math.floor((z-origin[2])/spacing); if(ix<0||ix>=dims[0]||iy<0||iy>=dims[1]||iz<0||iz>=dims[2]) return 0; return density[iz*dims[1]*dims[0]+iy*dims[0]+ix]; }}
  const resolution=80; const effect=new MarchingCubes(resolution,null,true,true,500000); effect.enableUvs=false; effect.enableColors=true;
  const isoMaterial=new THREE.MeshPhongMaterial({{color:0x60a5fa,transparent:true,opacity:0.75,shininess:80,specular:0x404080,side:THREE.DoubleSide,vertexColors:true}}); effect.material=isoMaterial;
  const size={{x:(dims[0]-1)*spacing,y:(dims[1]-1)*spacing,z:(dims[2]-1)*spacing}};
  effect.scale.set(size.x,size.y,size.z); effect.position.set(origin[0]+size.x/2,origin[1]+size.y/2,origin[2]+size.z/2);
  function updateIsosurface(level){{
    effect.reset(); const field=new Float32Array(resolution**3);
    for(let iz=0;iz<resolution;iz++){{ const z=origin[2]+(iz/(resolution-1))*size.z; for(let iy=0;iy<resolution;iy++){{ const y=origin[1]+(iy/(resolution-1))*size.y; for(let ix=0;ix<resolution;ix++){{ const x=origin[0]+(ix/(resolution-1))*size.x; field[iz*resolution*resolution+iy*resolution+ix]=getDensity(x,y,z)-level; }} }} }}
    effect.generate(field,level);
    const geometry=effect.geometry;
    if(geometry.attributes.position&&geometry.attributes.position.count>0){{ const colors=new Float32Array(geometry.attributes.position.count*3); for(let i=0;i<geometry.attributes.position.count;i++){{ colors[i*3]=0.376+0.2*Math.random(); colors[i*3+1]=0.647+0.15*Math.random(); colors[i*3+2]=0.980+0.02*Math.random(); }} geometry.setAttribute('color',new THREE.BufferAttribute(colors,3)); }}
  }}
  const moleculeGroup=new THREE.Group();
  MOLECULE_DATA.atoms.forEach(atom=>{{ const color=ELEMENT_COLORS[atom.element]||0xcccccc; const radius=Math.max(0.2,atom.radius*0.5); const geom=new THREE.SphereGeometry(radius,32,32); const mat=new THREE.MeshPhongMaterial({{color,shininess:120,specular:0x888888}}); const sphere=new THREE.Mesh(geom,mat); sphere.position.set(atom.x,atom.y,atom.z); sphere.castShadow=true; moleculeGroup.add(sphere); }});
  for(let i=0;i<MOLECULE_DATA.atoms.length;i++){{ for(let j=i+1;j<MOLECULE_DATA.atoms.length;j++){{ const a=MOLECULE_DATA.atoms[i],b=MOLECULE_DATA.atoms[j]; const dx=b.x-a.x,dy=b.y-a.y,dz=b.z-a.z; const dist=Math.sqrt(dx*dx+dy*dy+dz*dz); const cov=a.radius+b.radius; if(dist<cov*1.2){{ const dir=new THREE.Vector3(dx,dy,dz); const bondGeom=new THREE.CylinderGeometry(0.06,0.06,dist,12); const bondMat=new THREE.MeshPhongMaterial({{color:0x888899,shininess:100}}); const bond=new THREE.Mesh(bondGeom,bondMat); bond.position.set((a.x+b.x)/2,(a.y+b.y)/2,(a.z+b.z)/2); const q=new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0,1,0),dir.clone().normalize()); bond.quaternion.copy(q); moleculeGroup.add(bond); }} }} }}
  scene.add(moleculeGroup); scene.add(effect);
  const box=new THREE.Box3().setFromObject(moleculeGroup); const center=box.getCenter(new THREE.Vector3()); const radius=box.getSize(new THREE.Vector3()).length()*0.7; controls.target.copy(center); camera.position.copy(center).add(new THREE.Vector3(radius,radius*0.7,radius*1.2)); controls.update();
  let currentLevel=0.05; updateIsosurface(currentLevel);
  document.getElementById('isoLevel').addEventListener('input',e=>{{ currentLevel=parseFloat(e.target.value); document.getElementById('isoLevelValue').textContent=currentLevel.toFixed(3); updateIsosurface(currentLevel); }});
  document.getElementById('opacity').addEventListener('input',e=>{{ isoMaterial.opacity=parseFloat(e.target.value); document.getElementById('opacityValue').textContent=isoMaterial.opacity.toFixed(2); }});
  document.getElementById('showDensity').addEventListener('change',()=>{{ effect.visible=document.getElementById('showDensity').checked; }});
  document.getElementById('showAtoms').addEventListener('change',()=>{{ moleculeGroup.children.forEach(c=>{{ if(c.geometry instanceof THREE.SphereGeometry) c.visible=document.getElementById('showAtoms').checked; }}); }});
  document.getElementById('showBonds').addEventListener('change',()=>{{ moleculeGroup.children.forEach(c=>{{ if(c.geometry instanceof THREE.CylinderGeometry) c.visible=document.getElementById('showBonds').checked; }}); }});
  window.addEventListener('resize',()=>{{ camera.aspect=container.clientWidth/container.clientHeight; camera.updateProjectionMatrix(); renderer.setSize(container.clientWidth,container.clientHeight); }});
  (function animate(){{ requestAnimationFrame(animate); controls.update(); renderer.render(scene,camera); }})();
}}
initEnergyDiagram(); initInfoPanel(); initThreeJS();
</script>
</body></html>"#
    )
}

pub fn get_reaction_path_html(
    rp_json: &str,
) -> String {
    format!(
        r##"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>Reaction Path - Energy Evolution</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; background:linear-gradient(135deg,#0f0c29,#302b63,#24243e); color:#e0e0e0; min-height:100vh; overflow-x:hidden; }}
.header {{ padding:24px 32px; background:rgba(0,0,0,0.3); border-bottom:1px solid rgba(255,255,255,0.1); }}
.header h1 {{ font-size:24px; background:linear-gradient(90deg,#60a5fa,#a78bfa,#f472b6); -webkit-background-clip:text; -webkit-text-fill-color:transparent; margin-bottom:8px; }}
.header p {{ font-size:13px; color:#9ca3af; }}
.container {{ display:grid; grid-template-columns:1fr 1fr; gap:20px; padding:20px; }}
@media (max-width:1200px) {{ .container {{ grid-template-columns:1fr; }} }}
.panel {{ background:rgba(255,255,255,0.04); border:1px solid rgba(255,255,255,0.08); border-radius:12px; padding:20px; backdrop-filter:blur(10px); }}
.panel-title {{ font-size:16px; font-weight:600; color:#c4b5fd; margin-bottom:16px; padding-bottom:10px; border-bottom:1px solid rgba(196,181,253,0.2); display:flex; align-items:center; gap:8px; }}
.panel-title::before {{ content:''; width:4px; height:16px; background:linear-gradient(180deg,#60a5fa,#a78bfa); border-radius:2px; }}
#viewer3d {{ width:100%; height:480px; border-radius:8px; background:radial-gradient(ellipse at center,#1a1a2e 0%,#0a0a14 100%); }}
#energyChart {{ width:100%; height:480px; border-radius:8px; background:#111127; position:relative; }}
#energyCanvas {{ width:100%; height:100%; }}
.timeline-bar {{ background:rgba(167,139,250,0.06); border:1px solid rgba(167,139,250,0.12); border-radius:12px; padding:20px; margin:20px; }}
.timeline-controls {{ display:flex; align-items:center; gap:16px; flex-wrap:wrap; }}
.timeline-controls label {{ font-size:13px; color:#9ca3af; white-space:nowrap; }}
.timeline-controls input[type=range] {{ flex:1; min-width:200px; accent-color:#a78bfa; height:6px; }}
.step-display {{ font-size:22px; font-weight:700; color:#e0e7ff; font-family:'SF Mono',monospace; min-width:80px; text-align:center; }}
.play-btn {{ background:linear-gradient(135deg,#7c3aed,#6366f1); border:none; color:white; border-radius:8px; padding:8px 20px; font-size:14px; cursor:pointer; font-weight:600; }}
.play-btn:hover {{ filter:brightness(1.2); }}
.speed-select {{ background:#1e1b4b; color:#c4b5fd; border:1px solid rgba(167,139,250,0.3); border-radius:6px; padding:4px 8px; font-size:12px; }}
.info-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; margin-top:16px; }}
.info-card {{ background:rgba(167,139,250,0.08); border:1px solid rgba(167,139,250,0.15); border-radius:8px; padding:12px; }}
.info-label {{ font-size:11px; color:#9ca3af; text-transform:uppercase; letter-spacing:0.5px; }}
.info-value {{ font-size:16px; font-weight:700; color:#e0e7ff; margin-top:4px; font-family:'SF Mono',monospace; }}
.controls {{ display:flex; flex-wrap:wrap; gap:10px; margin-top:16px; }}
.control-group {{ display:flex; flex-direction:column; gap:6px; flex:1; min-width:120px; }}
.control-group label {{ font-size:12px; color:#9ca3af; }}
.control-group input[type=range] {{ width:100%; accent-color:#a78bfa; }}
.control-group .value-display {{ font-size:11px; color:#c4b5fd; text-align:right; font-family:monospace; }}
</style>
</head>
<body>
<div class="header">
  <h1>Reaction Path Simulation</h1>
  <p>IRC interpolation: electron density evolution & energy level dynamics</p>
</div>

<div class="container">
  <div class="panel">
    <div class="panel-title">3D Electron Density (Animated)</div>
    <div id="viewer3d"></div>
    <div class="controls">
      <div class="control-group"><label>Isosurface Level</label><input type="range" id="isoLevel" min="0.005" max="0.3" step="0.005" value="0.05"/><div class="value-display" id="isoLevelValue">0.050</div></div>
      <div class="control-group"><label>Opacity</label><input type="range" id="opacity" min="0.1" max="1.0" step="0.05" value="0.70"/><div class="value-display" id="opacityValue">0.70</div></div>
    </div>
    <div class="info-grid" id="infoGrid"></div>
  </div>
  <div class="panel">
    <div class="panel-title">Energy Level Evolution</div>
    <div id="energyChart"><canvas id="energyCanvas"></canvas></div>
  </div>
</div>

<div class="timeline-bar">
  <div class="timeline-controls">
    <button class="play-btn" id="playBtn">&#9654; Play</button>
    <label>Reaction Coordinate:</label>
    <input type="range" id="timeSlider" min="0" max="100" step="1" value="0"/>
    <div class="step-display" id="stepDisplay">t=0.00</div>
    <label>Speed:</label>
    <select class="speed-select" id="speedSelect">
      <option value="2000">0.5x</option>
      <option value="1000" selected>1x</option>
      <option value="500">2x</option>
      <option value="250">4x</option>
    </select>
  </div>
</div>

<script type="importmap">{{ "imports": {{ "three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/" }} }}</script>
<script type="module">
import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';
import {{ MarchingCubes }} from 'three/addons/objects/MarchingCubes.js';

const RP_DATA = {rp_json};

const ELEMENT_COLORS = {{ 'H':0xffffff,'C':0x303030,'N':0x3050f8,'O':0xff0d0d,'F':0x90e050,'P':0xff8000,'S':0xffff30,'Cl':0x1ff01f,'Pt':0xd0d0d0,'Au':0xffd123,'Fe':0xe06633,'Cu':0xb87333 }};

let currentStep = 0;
const totalSteps = RP_DATA.steps.length;
let playing = false;
let playTimer = null;

function b64ToF32(b64, count) {{
  const binary = atob(b64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return new Float32Array(bytes.buffer, 0, count);
}}

function initEnergyCanvas() {{
  const canvas = document.getElementById('energyCanvas');
  const container = document.getElementById('energyChart');
  canvas.width = container.clientWidth * 2;
  canvas.height = container.clientHeight * 2;
  canvas.style.width = container.clientWidth + 'px';
  canvas.style.height = container.clientHeight + 'px';
  drawEnergyChart(currentStep);
}}

function drawEnergyChart(highlightStep) {{
  const canvas = document.getElementById('energyCanvas');
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const steps = RP_DATA.steps;
  const maxEv = Math.max(...steps.flatMap(s => s.eigenvalues));
  const minEv = Math.min(...steps.flatMap(s => s.eigenvalues));
  const pad = (maxEv - minEv) * 0.1 || 1;
  const yMin = minEv - pad, yMax = maxEv + pad;
  const marginL = 80, marginR = 30, marginT = 30, marginB = 50;
  const plotW = W - marginL - marginR, plotH = H - marginT - marginB;

  ctx.fillStyle = '#111127';
  ctx.fillRect(0, 0, W, H);

  ctx.strokeStyle = 'rgba(255,255,255,0.08)';
  ctx.lineWidth = 1;
  for (let i = 0; i <= 10; i++) {{
    const y = marginT + plotH * i / 10;
    ctx.beginPath(); ctx.moveTo(marginL, y); ctx.lineTo(marginL + plotW, y); ctx.stroke();
  }}

  ctx.fillStyle = '#9ca3af'; ctx.font = '20px monospace'; ctx.textAlign = 'right';
  for (let i = 0; i <= 10; i++) {{
    const val = yMax - (yMax - yMin) * i / 10;
    const y = marginT + plotH * i / 10;
    ctx.fillText(val.toFixed(1), marginL - 8, y + 5);
  }}

  const numLevels = Math.min(steps[0]?.eigenvalues.length || 0, 30);
  const levelColors = [];
  for (let i = 0; i < numLevels; i++) {{
    const hue = (i / numLevels) * 270;
    levelColors.push('hsl(' + hue + ',80%,65%)');
  }}

  for (let lvl = 0; lvl < numLevels; lvl++) {{
    ctx.strokeStyle = levelColors[lvl];
    ctx.lineWidth = 2;
    ctx.globalAlpha = 0.6;
    ctx.beginPath();
    for (let s = 0; s < steps.length; s++) {{
      const x = marginL + (s / (steps.length - 1)) * plotW;
      const ev = steps[s].eigenvalues[lvl] || 0;
      const y = marginT + (1 - (ev - yMin) / (yMax - yMin)) * plotH;
      if (s === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();
  }}
  ctx.globalAlpha = 1.0;

  const nElec = RP_DATA.num_electrons;
  const homoIdx = Math.floor((nElec - 1) / 2);
  if (homoIdx < numLevels) {{
    ctx.strokeStyle = '#60a5fa'; ctx.lineWidth = 3; ctx.setLineDash([]);
    ctx.beginPath();
    for (let s = 0; s < steps.length; s++) {{
      const x = marginL + (s / (steps.length - 1)) * plotW;
      const ev = steps[s].eigenvalues[homoIdx] || 0;
      const y = marginT + (1 - (ev - yMin) / (yMax - yMin)) * plotH;
      if (s === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();
  }}
  if (homoIdx + 1 < numLevels) {{
    ctx.strokeStyle = '#f472b6'; ctx.lineWidth = 3; ctx.setLineDash([8, 4]);
    ctx.beginPath();
    for (let s = 0; s < steps.length; s++) {{
      const x = marginL + (s / (steps.length - 1)) * plotW;
      const ev = steps[s].eigenvalues[homoIdx + 1] || 0;
      const y = marginT + (1 - (ev - yMin) / (yMax - yMin)) * plotH;
      if (s === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();
    ctx.setLineDash([]);
  }}

  const hlX = marginL + (highlightStep / (steps.length - 1)) * plotW;
  ctx.strokeStyle = 'rgba(250,204,21,0.5)'; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(hlX, marginT); ctx.lineTo(hlX, marginT + plotH); ctx.stroke();

  ctx.fillStyle = 'rgba(250,204,21,0.15)';
  ctx.fillRect(hlX - 6, marginT, 12, plotH);

  for (let lvl = 0; lvl < Math.min(numLevels, 30); lvl++) {{
    const ev = steps[highlightStep]?.eigenvalues[lvl] || 0;
    const y = marginT + (1 - (ev - yMin) / (yMax - yMin)) * plotH;
    ctx.fillStyle = levelColors[lvl]; ctx.beginPath(); ctx.arc(hlX, y, 5, 0, Math.PI * 2); ctx.fill();
  }}

  ctx.fillStyle = '#9ca3af'; ctx.font = '20px monospace'; ctx.textAlign = 'center';
  ctx.fillText('Reaction Coordinate (t)', marginL + plotW / 2, H - 8);

  const totalE = steps.map(s => s.total_energy);
  const minE = Math.min(...totalE), maxE = Math.max(...totalE);
  const ePad = (maxE - minE) * 0.1 || 1;
  ctx.strokeStyle = '#fbbf24'; ctx.lineWidth = 2;
  ctx.beginPath();
  for (let s = 0; s < steps.length; s++) {{
    const x = marginL + (s / (steps.length - 1)) * plotW;
    const norm = (totalE[s] - minE + ePad) / (maxE - minE + 2 * ePad);
    const y = marginT + (1 - norm) * plotH;
    if (s === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }}
  ctx.stroke();
}}

let scene, camera, renderer, controls, effect, isoMaterial, moleculeGroup;

function initThreeJS() {{
  const container = document.getElementById('viewer3d');
  scene = new THREE.Scene(); scene.background = new THREE.Color(0x0a0a14); scene.fog = new THREE.FogExp2(0x0a0a14, 0.015);
  camera = new THREE.PerspectiveCamera(60, container.clientWidth / container.clientHeight, 0.1, 1000); camera.position.set(8, 6, 12);
  renderer = new THREE.WebGLRenderer({{ antialias: true }}); renderer.setSize(container.clientWidth, container.clientHeight); renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(renderer.domElement);
  controls = new OrbitControls(camera, renderer.domElement); controls.enableDamping = true; controls.dampingFactor = 0.08;
  scene.add(new THREE.AmbientLight(0x404050, 0.6));
  const keyLight = new THREE.DirectionalLight(0xffffff, 1.2); keyLight.position.set(10, 15, 8); scene.add(keyLight);
  const fillLight = new THREE.DirectionalLight(0x8080ff, 0.5); fillLight.position.set(-8, 5, -6); scene.add(fillLight);
  scene.add(new THREE.PointLight(0xff80ff, 0.6, 50).position.set(-5, -3, -8) || new THREE.PointLight(0xff80ff, 0.6, 50));

  isoMaterial = new THREE.MeshPhongMaterial({{ color: 0x60a5fa, transparent: true, opacity: 0.7, shininess: 80, specular: 0x404080, side: THREE.DoubleSide, vertexColors: true }});
  effect = new MarchingCubes(64, null, true, true, 300000);
  effect.enableUvs = false; effect.enableColors = true; effect.material = isoMaterial;

  moleculeGroup = new THREE.Group();
  scene.add(moleculeGroup); scene.add(effect);

  updateStep(0);

  document.getElementById('isoLevel').addEventListener('input', e => {{
    document.getElementById('isoLevelValue').textContent = parseFloat(e.target.value).toFixed(3);
    updateStep(currentStep);
  }});
  document.getElementById('opacity').addEventListener('input', e => {{
    isoMaterial.opacity = parseFloat(e.target.value);
    document.getElementById('opacityValue').textContent = isoMaterial.opacity.toFixed(2);
  }});

  window.addEventListener('resize', () => {{
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  }});

  (function animate() {{ requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }})();
}}

function updateStep(stepIdx) {{
  const step = RP_DATA.steps[stepIdx];
  if (!step) return;
  currentStep = stepIdx;

  const dims = RP_DATA.grid_dims;
  const origin = RP_DATA.origin;
  const spacing = RP_DATA.spacing;
  const size = {{ x: (dims[0]-1)*spacing, y: (dims[1]-1)*spacing, z: (dims[2]-1)*spacing }};

  const density = b64ToF32(step.density_b64, dims[0]*dims[1]*dims[2]);

  const resolution = 64;
  effect.reset();
  effect.scale.set(size.x, size.y, size.z);
  effect.position.set(origin[0]+size.x/2, origin[1]+size.y/2, origin[2]+size.z/2);

  const isoLevel = parseFloat(document.getElementById('isoLevel').value);
  const field = new Float32Array(resolution ** 3);
  for (let iz = 0; iz < resolution; iz++) {{
    const z = origin[2] + (iz/(resolution-1))*size.z;
    for (let iy = 0; iy < resolution; iy++) {{
      const y = origin[1] + (iy/(resolution-1))*size.y;
      for (let ix = 0; ix < resolution; ix++) {{
        const x = origin[0] + (ix/(resolution-1))*size.x;
        const dIx = Math.floor((x-origin[0])/spacing);
        const dIy = Math.floor((y-origin[1])/spacing);
        const dIz = Math.floor((z-origin[2])/spacing);
        let d = 0;
        if (dIx >= 0 && dIx < dims[0] && dIy >= 0 && dIy < dims[1] && dIz >= 0 && dIz < dims[2]) {{
          d = density[dIz*dims[1]*dims[0] + dIy*dims[0] + dIx];
        }}
        field[iz*resolution*resolution + iy*resolution + ix] = d - isoLevel;
      }}
    }}
  }}
  effect.generate(field, isoLevel);

  while (moleculeGroup.children.length > 0) moleculeGroup.remove(moleculeGroup.children[0]);
  step.atoms.forEach(atom => {{
    const color = ELEMENT_COLORS[atom.element] || 0xcccccc;
    const radius = Math.max(0.2, (atom.radius || 0.7) * 0.5);
    const geom = new THREE.SphereGeometry(radius, 24, 24);
    const mat = new THREE.MeshPhongMaterial({{ color, shininess: 120 }});
    const mesh = new THREE.Mesh(geom, mat);
    mesh.position.set(atom.x, atom.y, atom.z);
    moleculeGroup.add(mesh);
  }});

  if (step.bonds) {{
    step.bonds.forEach(bond => {{
      const a = step.atoms[bond.atom1], b2 = step.atoms[bond.atom2];
      if (!a || !b2) return;
      const dx = b2.x-a.x, dy = b2.y-a.y, dz = b2.z-a.z;
      const dist = Math.sqrt(dx*dx+dy*dy+dz*dz);
      if (dist < 0.1) return;
      const dir = new THREE.Vector3(dx, dy, dz);
      const bGeom = new THREE.CylinderGeometry(0.04*bond.bond_order, 0.04*bond.bond_order, dist, 8);
      const bMat = new THREE.MeshPhongMaterial({{ color: 0x888899 }});
      const bMesh = new THREE.Mesh(bGeom, bMat);
      bMesh.position.set((a.x+b2.x)/2,(a.y+b2.y)/2,(a.z+b2.z)/2);
      bMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), dir.clone().normalize());
      moleculeGroup.add(bMesh);
    }});
  }}

  document.getElementById('stepDisplay').textContent = 't=' + step.t.toFixed(2);
  document.getElementById('timeSlider').value = stepIdx;

  const grid = document.getElementById('infoGrid');
  const nElec = RP_DATA.num_electrons;
  const homoI = Math.floor((nElec-1)/2);
  const cards = [
    {{ label:'Step', value: step.step + '/' + (totalSteps-1) }},
    {{ label:'t', value: step.t.toFixed(3) }},
    {{ label:'Total Energy', value: step.total_energy.toFixed(2) + ' eV' }},
    {{ label:'HOMO-LUMO Gap', value: step.homo_lumo_gap.toFixed(4) + ' eV' }},
    {{ label:'HOMO (eV)', value: (step.eigenvalues[homoI]||0).toFixed(4) }},
    {{ label:'Orbitals', value: step.eigenvalues.length }},
  ];
  grid.innerHTML = cards.map(c => '<div class="info-card"><div class="info-label">'+c.label+'</div><div class="info-value">'+c.value+'</div></div>').join('');

  drawEnergyChart(stepIdx);
}}

document.getElementById('timeSlider').addEventListener('input', e => {{
  const step = parseInt(e.target.value);
  if (step >= 0 && step < totalSteps) updateStep(step);
}});

document.getElementById('playBtn').addEventListener('click', () => {{
  playing = !playing;
  const btn = document.getElementById('playBtn');
  btn.innerHTML = playing ? '&#9646;&#9646; Pause' : '&#9654; Play';
  if (playing) {{
    if (currentStep >= totalSteps - 1) currentStep = 0;
    startPlayback();
  }} else {{
    stopPlayback();
  }}
}});

function startPlayback() {{
  stopPlayback();
  const speed = parseInt(document.getElementById('speedSelect').value);
  playTimer = setInterval(() => {{
    if (currentStep < totalSteps - 1) {{
      updateStep(currentStep + 1);
    }} else {{
      playing = false;
      document.getElementById('playBtn').innerHTML = '&#9654; Play';
      stopPlayback();
    }}
  }}, speed);
}}

function stopPlayback() {{
  if (playTimer) {{ clearInterval(playTimer); playTimer = null; }}
}}

initThreeJS();
initEnergyCanvas();

window.addEventListener('resize', () => {{
  const canvas = document.getElementById('energyCanvas');
  const container = document.getElementById('energyChart');
  canvas.width = container.clientWidth * 2;
  canvas.height = container.clientHeight * 2;
  canvas.style.width = container.clientWidth + 'px';
  canvas.style.height = container.clientHeight + 'px';
  drawEnergyChart(currentStep);
}});
</script>
</body></html>"##
    )
}
