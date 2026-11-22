/* ===== Laboratorio 3D de simetría (Three.js r160 + WebGL2) =====
   Depende de SYM (sym.js) y ORB_SHADERS (shaders.js), cargados como scripts clásicos antes. */
const host = document.getElementById('s3View');
const msg = t => { host.innerHTML = `<div class="s3msg">${t}</div>`; };
let THREE, OrbitControls;
try {
  THREE = await import('three');
  ({ OrbitControls } = await import('three/addons/controls/OrbitControls.js'));
} catch (e) {
  try { sym3dLite(host); } catch (e2) { msg('Sin conexión: no se pudo cargar Three.js y la vista simplificada falló.'); }
}
if (THREE) init();

function init(){
  const $ = s => document.querySelector(s);
  const f3 = x => (Math.abs(x) < 5e-4 ? 0 : x).toFixed(3).replace('-', '−');
  const canvas = document.createElement('canvas'); host.appendChild(canvas);

  /* ---------- renderer, cámara, controles ---------- */
  let renderer;
  try { renderer = new THREE.WebGLRenderer({ canvas, antialias: true, powerPreference: 'high-performance' }); }
  catch (e) { msg('Tu navegador no pudo crear un contexto WebGL2.'); return; }
  if (!renderer.capabilities.isWebGL2) { msg('Este laboratorio necesita WebGL2.'); return; }
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.5));
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0b1018);
  const camera = new THREE.PerspectiveCamera(38, 1, 0.1, 100);
  camera.up.set(0, 0, 1);                       // convención química: z hacia arriba
  camera.position.set(8.5, -10, 6.5);
  const controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true; controls.dampingFactor = 0.08;
  controls.minDistance = 4; controls.maxDistance = 32;
  scene.add(new THREE.AmbientLight(0xffffff, 0.55));
  const sun = new THREE.DirectionalLight(0xffffff, 1.4); sun.position.set(5, -4, 9); scene.add(sun);
  const fill = new THREE.DirectionalLight(0x8fb4ff, 0.5); fill.position.set(-6, 5, -3); scene.add(fill);

  /* ---------- estado ---------- */
  const S = { g: 'Oh', mode: 'd', idx: 2, overlay: true, acc: SYM.I3(), anim: null, speed: 1,
              showElems: true, showNodes: true, showPoly: true, showFaces: true, showLabels: true, log: [] };

  /* ---------- utilidades Three ---------- */
  const v3 = a => new THREE.Vector3(a[0], a[1], a[2]);
  const m3 = A => new THREE.Matrix3().set(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]);
  function label(text, color = '#e8eef7', size = 0.55){
    const c = document.createElement('canvas'); c.width = 128; c.height = 64;
    const x = c.getContext('2d'); x.font = '600 34px IBM Plex Sans, system-ui, sans-serif'; x.textAlign = 'center'; x.textBaseline = 'middle';
    x.lineWidth = 6; x.strokeStyle = 'rgba(8,12,20,.85)'; x.strokeText(text, 64, 32); x.fillStyle = color; x.fillText(text, 64, 32);
    const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false, transparent: true }));
    sp.scale.set(size * 2, size, 1); sp.renderOrder = 30; return sp;
  }
  function orient(obj, from, to){ obj.quaternion.setFromUnitVectors(from.clone().normalize(), to.clone().normalize()); }
  function disposeGroup(g){ g.traverse(o => { if (o.geometry) o.geometry.dispose(); if (o.material) { (Array.isArray(o.material) ? o.material : [o.material]).forEach(m => { if (m.map) m.map.dispose(); m.dispose(); }); } }); g.clear(); }

  /* ---------- ejes cartesianos ---------- */
  const axes = new THREE.Group(); scene.add(axes);
  [['x', [1,0,0], 0xd45d5d], ['y', [0,1,0], 0x5dbb73], ['z', [0,0,1], 0x5d8fd4]].forEach(([n, d, c]) => {
    const g = new THREE.BufferGeometry().setFromPoints([v3(d).multiplyScalar(-4.6), v3(d).multiplyScalar(4.6)]);
    axes.add(new THREE.Line(g, new THREE.LineBasicMaterial({ color: c, transparent: true, opacity: 0.55 })));
    const l = label(n, '#' + c.toString(16).padStart(6, '0'), 0.45); l.position.copy(v3(d).multiplyScalar(4.95)); axes.add(l);
  });

  /* ---------- complejo (se transforma con la operación) ---------- */
  const mol = new THREE.Group(); mol.matrixAutoUpdate = false; scene.add(mol);
  const nodes = new THREE.Group(); nodes.matrixAutoUpdate = false; scene.add(nodes);
  const elems = new THREE.Group(); scene.add(elems);
  const active = new THREE.Group(); scene.add(active);

  function buildMolecule(){
    disposeGroup(mol);
    const G = SYM.GROUPS[S.g], D = SYM.D;
    const metal = new THREE.Mesh(new THREE.SphereGeometry(0.26, 32, 16), new THREE.MeshStandardMaterial({ color: 0xb9c3d1, metalness: 0.6, roughness: 0.3 }));
    mol.add(metal);
    const ligMat = new THREE.MeshStandardMaterial({ color: 0x9fb0c8, roughness: 0.45 });
    const ligMat1 = new THREE.MeshStandardMaterial({ color: 0xf2b84b, roughness: 0.4, emissive: 0x3a2a00 });
    const bondMat = new THREE.MeshStandardMaterial({ color: 0x5f6f86, roughness: 0.6 });
    G.ligs.forEach((L, k) => {
      const p = v3(L).multiplyScalar(D);
      const s = new THREE.Mesh(new THREE.SphereGeometry(0.3, 28, 14), k === 0 ? ligMat1 : ligMat); s.position.copy(p); mol.add(s);
      const b = new THREE.Mesh(new THREE.CylinderGeometry(0.05, 0.05, D - 0.3, 10), bondMat);
      b.position.copy(p.clone().multiplyScalar(0.5 - 0.15 / D)); orient(b, new THREE.Vector3(0, 1, 0), p); mol.add(b);
      if (S.showLabels) { const l = label('L' + (k + 1), k === 0 ? '#f2b84b' : '#dfe7f2', 0.42); l.position.copy(p.clone().multiplyScalar(1.22)); mol.add(l); }
    });
    if (S.showPoly) {
      const pts = G.ligs.map(L => v3(L).multiplyScalar(D));
      const seg = []; G.edges.forEach(([a, b]) => seg.push(pts[a], pts[b]));
      mol.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(seg), new THREE.LineBasicMaterial({ color: 0x8aa2c4, transparent: true, opacity: 0.7 })));
    }
    if (S.showFaces) {
      // Caras translúcidas, un color por cara; depthWrite = false para que los átomos de detrás sigan viéndose
      let face;
      if (S.g === 'Oh') face = new THREE.OctahedronGeometry(D);
      else if (S.g === 'Td') face = new THREE.TetrahedronGeometry(D);           // vértices (±1,±1,±1) con producto +1: coinciden con los ligandos
      else { face = new THREE.PlaneGeometry(D * Math.SQRT2, D * Math.SQRT2); face.rotateZ(Math.PI / 4); face = face.toNonIndexed(); }
      const pal = [0x0072b2, 0xd55e00, 0x009e73, 0xe69f00, 0x56b4e9, 0xcc79a7], cnt = face.attributes.position.count, cols = new Float32Array(cnt * 3), tmp = new THREE.Color();
      for (let v = 0; v < cnt; v++) { tmp.setHex(pal[Math.floor(v / (S.g === 'D4h' ? 6 : 3)) % pal.length]); cols[3*v] = tmp.r; cols[3*v+1] = tmp.g; cols[3*v+2] = tmp.b; }
      face.setAttribute('color', new THREE.Float32BufferAttribute(cols, 3));
      const fm = new THREE.Mesh(face, new THREE.MeshStandardMaterial({ vertexColors: true, transparent: true, opacity: 0.24, side: THREE.DoubleSide, depthWrite: false,
        roughness: 0.55, metalness: 0, flatShading: true, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1 }));
      fm.renderOrder = 5; mol.add(fm);
    }
  }

  /* ---------- planos nodales del orbital base ---------- */
  function buildNodes(){
    disposeGroup(nodes);
    if (S.mode !== 'd') return;
    const mat = new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.09, side: THREE.DoubleSide, depthWrite: false });
    const edge = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.35 });
    const plane = n => { const g = new THREE.PlaneGeometry(7.5, 7.5); const m = new THREE.Mesh(g, mat); orient(m, new THREE.Vector3(0, 0, 1), v3(n)); m.renderOrder = 20;
      const e = new THREE.LineSegments(new THREE.EdgesGeometry(g), edge); e.quaternion.copy(m.quaternion); nodes.add(m, e); };
    const P = { 1: [[1,-1,0],[1,1,0]], 2: [[1,0,0],[0,1,0]], 3: [[1,0,0],[0,0,1]], 4: [[0,1,0],[0,0,1]] };
    if (S.idx === 0) {                                     // d_z²: superficie nodal cónica, θ = 54,74°
      const h = 3.6, rad = h * Math.SQRT2;
      [1, -1].forEach(sg => { const g = new THREE.ConeGeometry(rad, h, 64, 1, true); g.translate(0, -h / 2, 0); g.rotateX(sg > 0 ? -Math.PI / 2 : Math.PI / 2);
        const m = new THREE.Mesh(g, mat); m.renderOrder = 20; nodes.add(m); });
    } else P[S.idx].forEach(plane);
  }

  /* ---------- elementos de simetría ---------- */
  const ECOL = { 4: 0x5aa9ff, 3: 0x59d18c, 2: 0xf2b84b, h: 0x4dd0e1, v: 0x9ccc65, d: 0xf48fb1, i: 0xffffff };
  function axisMesh(n, k, strong, prime){
    const g = new THREE.Group(), col = ECOL[k] || ECOL[k / 2] || 0xffd36b;
    const c = new THREE.Mesh(new THREE.CylinderGeometry(strong ? 0.045 : 0.022, strong ? 0.045 : 0.022, 9.6, 8),
      new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: strong ? 0.95 : 0.55 }));
    orient(c, new THREE.Vector3(0, 1, 0), v3(n)); g.add(c);
    const l = label('C' + k + (prime === 1 ? '′' : prime === 2 ? '″' : ''), '#' + col.toString(16).padStart(6, '0'), strong ? 0.5 : 0.38);
    l.position.copy(v3(n).normalize().multiplyScalar(5.05)); g.add(l); return g;
  }
  function planeMesh(n, kind, strong){
    const g = new THREE.Group(), col = ECOL[kind];
    const geo = new THREE.PlaneGeometry(8.4, 8.4);
    const m = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: strong ? 0.2 : 0.045, side: THREE.DoubleSide, depthWrite: false }));
    orient(m, new THREE.Vector3(0, 0, 1), v3(n)); m.renderOrder = 6; g.add(m);
    const e = new THREE.LineSegments(new THREE.EdgesGeometry(geo), new THREE.LineBasicMaterial({ color: col, transparent: true, opacity: strong ? 0.9 : 0.25 }));
    e.quaternion.copy(m.quaternion); g.add(e); return g;
  }
  function invMesh(strong){ return new THREE.Mesh(new THREE.SphereGeometry(strong ? 0.42 : 0.34, 16, 10), new THREE.MeshBasicMaterial({ color: ECOL.i, wireframe: true, transparent: true, opacity: strong ? 0.9 : 0.4 })); }
  /* catálogo completo de elementos de simetría, con las operaciones que genera cada uno */
  const DN = n => { if (n.filter(x => x).length === 1) return '(' + 'xyz'[n.findIndex(x => x)] + ')'; return '[' + n.map(x => x < 0 ? '−1' : String(x)).join('') + ']'; };
  function catalog(g){
    const A = [], add = (grp, t, n, extra, ops) => A.push(Object.assign({ grp, t, n }, extra, { ops }));
    const ax = [[1,0,0],[0,1,0],[0,0,1]], c3 = [[1,1,1],[1,1,-1],[1,-1,1],[-1,1,1]], c2p = [[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]];
    if (g === 'Oh') {
      ax.forEach(n => add('Ejes C₄ (cada uno contiene también C₂ = C₄² y S₄)', 'C', n, { k: 4 }, [['C₄', 'C', 4, 3], ['C₂', 'C', 2, 4], ['S₄', 'S', 4, 6]]));
      c3.forEach(n => add('Ejes C₃ (cada uno contiene también S₆)', 'C', n, { k: 3 }, [['C₃', 'C', 3, 1], ['S₆', 'S', 6, 7]]));
      c2p.forEach(n => add('Ejes C₂′ (bisectan aristas opuestas)', 'C', n, { k: 2, p: 1 }, [['C₂′', 'C', 2, 2]]));
      ax.forEach(n => add('Planos σₕ (contienen 4 ligandos)', 'sigma', n, { kind: 'h' }, [['σₕ', 'sigma', 0, 8, 'sh']]));
      c2p.forEach(n => add('Planos σd (contienen 2 ligandos opuestos y bisectan)', 'sigma', n, { kind: 'd' }, [['σd', 'sigma', 0, 9, 'sd']]));
      add('Centro de inversión', 'i', null, {}, [['i', 'i', 0, 5]]);
    } else if (g === 'Td') {
      c3.forEach(n => add('Ejes C₃ (pasan por un ligando y el centro de la cara opuesta)', 'C', n, { k: 3 }, [['C₃', 'C', 3, 1]]));
      ax.forEach(n => add('Ejes C₂ (cada uno es también eje S₄)', 'C', n, { k: 2, s4: 1 }, [['C₂', 'C', 2, 2], ['S₄', 'S', 4, 3]]));
      c2p.forEach(n => add('Planos σd (contienen 2 ligandos)', 'sigma', n, { kind: 'd' }, [['σd', 'sigma', 0, 4, 'sd']]));
    } else {
      add('Eje principal C₄ (contiene también C₂ y S₄)', 'C', [0,0,1], { k: 4 }, [['C₄', 'C', 4, 1], ['C₂', 'C', 2, 2], ['S₄', 'S', 4, 6]]);
      [[1,0,0],[0,1,0]].forEach(n => add('Ejes C₂′ (pasan por los ligandos)', 'C', n, { k: 2, p: 1 }, [['C₂′', 'C', 2, 3]]));
      [[1,1,0],[1,-1,0]].forEach(n => add('Ejes C₂″ (entre ligandos)', 'C', n, { k: 2, p: 2 }, [['C₂″', 'C', 2, 4]]));
      add('Plano σₕ (plano molecular)', 'sigma', [0,0,1], { kind: 'h' }, [['σₕ', 'sigma', 0, 7, 'sh']]);
      [[0,1,0],[1,0,0]].forEach(n => add('Planos σᵥ (contienen ligandos y el eje C₄)', 'sigma', n, { kind: 'v' }, [['σᵥ', 'sigma', 0, 8, 'sv']]));
      [[1,1,0],[1,-1,0]].forEach(n => add('Planos σd (entre ligandos)', 'sigma', n, { kind: 'd' }, [['σd', 'sigma', 0, 9, 'sd']]));
      add('Centro de inversión', 'i', null, {}, [['i', 'i', 0, 5]]);
    }
    A.forEach(e => { e.name = e.t === 'i' ? 'i (centro)' : (e.t === 'C' ? 'C' + '₀₁₂₃₄₅₆'[e.k] + (e.p === 1 ? '′' : e.p === 2 ? '″' : '') : 'σ' + ({ h: 'ₕ', v: 'ᵥ', d: 'd' }[e.kind])) + (e.n ? ' ' + DN(e.n) : '');
      e.opObjs = e.ops.map(([lb, type, k, cls, id]) => ({ id: id || (type + k), label: lb + (e.n ? ' ' + DN(e.n) : ''), type, n: e.n, k, cls })); });
    return A;
  }
  const ECSS = k => '#' + ECOL[k].toString(16).padStart(6, '0');
  function buildElems(){
    disposeGroup(elems);
    if (!S.showElems) return;
    (S.cat || []).forEach((e, i) => {
      if (!S.eOn[i]) return; const strong = S.hover === i;
      if (e.t === 'C') { const a = axisMesh(e.n, e.k, strong, e.p); elems.add(a);
        if (e.s4 || e.k === 4) { const l = label('S4', '#ffd36b', 0.3); l.position.copy(v3(e.n).normalize().multiplyScalar(-5.05)); elems.add(l); }
        if (e.k === 3 && S.g === 'Oh') { const l = label('S6', '#ffd36b', 0.3); l.position.copy(v3(e.n).normalize().multiplyScalar(-5.05)); elems.add(l); } }
      else if (e.t === 'sigma') elems.add(planeMesh(e.n, e.kind, strong));
      else elems.add(invMesh(strong));
    });
  }
  function buildEList(){
    const box = $('#s3EList'); if (!box) return;
    const G = SYM.GROUPS[S.g], cat = S.cat, groups = [];
    cat.forEach((e, i) => { let g = groups.find(x => x.n === e.grp); if (!g) groups.push(g = { n: e.grp, idx: [] }); g.idx.push(i); });
    const sizes = G.classes.map(c => { const m = c.match(/^(\d+)/); return m ? +m[1] : 1; }), h = sizes.reduce((a, b) => a + b, 0);
    const col = e => e.t === 'C' ? ECSS(e.k) : e.t === 'sigma' ? ECSS(e.kind) : '#ffffff';
    box.innerHTML = `<div class="grp" style="margin-top:0"><span>${cat.length} elementos · orden del grupo h = ${h} operaciones</span><span><a data-all="1">todos</a><a data-all="0">ninguno</a></span></div>` +
      groups.map(g => `<div class="grp"><span>${g.n} · ${g.idx.length}</span><span><a data-only="${g.idx.join(',')}">solo estos</a></span></div>` +
        g.idx.map(i => { const e = cat[i]; return `<div class="row" data-i="${i}"><input type="checkbox" data-e="${i}" ${S.eOn[i] ? 'checked' : ''} aria-label="Mostrar ${e.name}"><i style="display:inline-block;width:10px;height:10px;border-radius:2px;background:${col(e)}"></i><span style="flex:1">${e.name}</span>` +
          e.opObjs.map((o, j) => `<button data-i="${i}" data-j="${j}" title="Aplicar ${o.label} · clase ${G.classes[o.cls]}">${e.ops[j][0]}</button>`).join('') + '</div>'; }).join('')).join('');
    box.querySelectorAll('input[data-e]').forEach(c => c.addEventListener('change', () => { S.eOn[+c.dataset.e] = c.checked; buildElems(); }));
    box.querySelectorAll('button[data-j]').forEach(b => b.addEventListener('click', () => { const e = cat[+b.dataset.i]; S.eOn[+b.dataset.i] = true; if (!S.showElems) { S.showElems = true; $('#s3Elems').checked = true; } buildEList(); buildElems(); runOp(e.opObjs[+b.dataset.j]); }));
    box.querySelectorAll('a[data-all]').forEach(a => a.addEventListener('click', () => { S.eOn = cat.map(() => a.dataset.all === '1'); buildEList(); buildElems(); }));
    box.querySelectorAll('a[data-only]').forEach(a => a.addEventListener('click', () => { const keep = a.dataset.only.split(',').map(Number); S.eOn = cat.map((_, i) => keep.includes(i)); buildEList(); buildElems(); }));
    box.querySelectorAll('.row').forEach(r => { r.addEventListener('mouseenter', () => { S.hover = +r.dataset.i; buildElems(); }); r.addEventListener('mouseleave', () => { S.hover = null; buildElems(); }); });
    box.style.opacity = S.showElems ? 1 : 0.45;
  }
  function buildActive(op){
    disposeGroup(active);
    if (!op || op.type === 'E') return;
    if (op.type === 'C') active.add(axisMesh(op.n, op.k, true));
    if (op.type === 'sigma') active.add(planeMesh(op.n, op.id === 'sh' ? 'h' : op.id === 'sv' ? 'v' : 'd', true));
    if (op.type === 'i') active.add(invMesh(true));
    if (op.type === 'S') { const a = axisMesh(op.n, op.k, true); a.children[1].material.map.dispose(); a.remove(a.children[1]);
      const l = label('S' + op.k, '#ffd36b', 0.5); l.position.copy(v3(op.n).normalize().multiplyScalar(5.05)); a.add(l);
      active.add(a, planeMesh(op.n, 'h', true)); }
  }

  /* ---------- orbital por raymarching ---------- */
  const BOX = 4.6;
  const uniforms = {
    uCam: { value: new THREE.Vector3() }, uViewProj: { value: new THREE.Matrix4() }, uMinv: { value: new THREE.Matrix3() },
    uC: { value: new Float32Array(17) }, uLig: { value: Array.from({ length: 8 }, () => new THREE.Vector3(0, 0, 99)) },
    uPiV: { value: Array.from({ length: 8 }, () => new THREE.Vector3()) },
    uIso: { value: 0.3 }, uOpacity: { value: 0.85 }, uColPos: { value: new THREE.Vector3(0.31, 0.525, 1.0) }, uColNeg: { value: new THREE.Vector3(1.0, 0.54, 0.24) },   // sRGB directo: el shader crudo no convierte
    uBox: { value: BOX }, uSteps: { value: 160 }
  };
  const orbMat = new THREE.RawShaderMaterial({ vertexShader: ORB_SHADERS.vert, fragmentShader: ORB_SHADERS.frag, glslVersion: THREE.GLSL3,
    uniforms, transparent: true, side: THREE.BackSide, depthWrite: true, depthTest: true });
  const orb = new THREE.Mesh(new THREE.BoxGeometry(2 * BOX, 2 * BOX, 2 * BOX), orbMat);
  orb.renderOrder = 10; orb.frustumCulled = false; scene.add(orb);
  function setFunction(){
    const f = SYM.renderVector(S.g, S.mode, S.idx, S.overlay);
    uniforms.uC.value.set(f.slice(0, 17));
    uniforms.uPiV.value.forEach((v, k) => v.set(f[17 + 3*k], f[18 + 3*k], f[19 + 3*k]));
    const L = SYM.GROUPS[S.g].ligs;
    uniforms.uLig.value.forEach((v, k) => k < L.length ? v.set(L[k][0] * SYM.D, L[k][1] * SYM.D, L[k][2] * SYM.D) : v.set(0, 0, 99));
  }

  /* ---------- interfaz ---------- */
  const seg = (sel, items, cur, cb) => { const el = $(sel); el.innerHTML = items.map(([v, t]) => `<button data-v="${v}" aria-pressed="${String(v) === String(cur)}">${t}</button>`).join('');
    el.querySelectorAll('button').forEach(b => b.addEventListener('click', () => { el.querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', 'false')); b.setAttribute('aria-pressed', 'true'); cb(b.dataset.v); })); };
  const ORBN = ['d<sub>z²</sub>', 'd<sub>x²−y²</sub>', 'd<sub>xy</sub>', 'd<sub>xz</sub>', 'd<sub>yz</sub>'];
  function buildUI(){
    const G = SYM.GROUPS[S.g];
    seg('#s3Group', [['Oh', 'O<sub>h</sub> · octaedro'], ['Td', 'T<sub>d</sub> · tetraedro'], ['D4h', 'D<sub>4h</sub> · plano cuadrado']], S.g, v => { S.g = v; S.idx = S.mode === 'd' ? S.idx : 0; resetAll(); });
    seg('#s3Mode', [['d', 'Orbital d del metal'], ['salc', 'SALC σ de ligandos']], S.mode, v => { S.mode = v; S.idx = v === 'd' ? 2 : 0; resetAll(); });
    if (S.mode === 'd') { seg('#s3Fn', ORBN.map((t, i) => [i, t]), S.idx, v => { S.idx = +v; resetAll(); }); $('#s3OverlayRow').hidden = true; }
    else { seg('#s3Fn', G.salcs.map((s, i) => [i, s.name.replace('_g', '<sub>g</sub>').replace('_u', '<sub>u</sub>')]), S.idx, v => { S.idx = +v; resetAll(); }); $('#s3OverlayRow').hidden = false; }
    const ops = $('#s3Ops'); ops.innerHTML = G.ops.map((o, i) => `<button data-i="${i}" title="Clase ${G.classes[o.cls]}">${o.label}</button>`).join('') + '<button id="s3Reset" class="reset">Reiniciar</button>';
    ops.querySelectorAll('button[data-i]').forEach(b => b.addEventListener('click', () => runOp(G.ops[+b.dataset.i])));
    if (S.catG !== S.g) { S.cat = catalog(S.g); S.eOn = S.cat.map(() => true); S.catG = S.g; S.hover = null; }
    buildEList();
    $('#s3Reset').addEventListener('click', resetAll);
  }
  function resetAll(){ S.acc = SYM.I3(); S.anim = null; S.log = []; buildUI(); buildMolecule(); buildNodes(); buildElems(); buildActive(null); setFunction(); describe(null); }
  $('#s3Overlay').addEventListener('change', e => { S.overlay = e.target.checked; setFunction(); });
  $('#s3Elems').addEventListener('change', e => { S.showElems = e.target.checked; buildElems(); buildEList(); });
  $('#s3Nodes').addEventListener('change', e => { S.showNodes = e.target.checked; nodes.visible = S.showNodes; });
  $('#s3Poly').addEventListener('change', e => { S.showPoly = e.target.checked; buildMolecule(); });
  $('#s3Faces').addEventListener('change', e => { S.showFaces = e.target.checked; buildMolecule(); });
  $('#s3Labels').addEventListener('change', e => { S.showLabels = e.target.checked; buildMolecule(); });
  const slider = (id, out, fn, fmt) => { const el = $(id); const upd = () => { const v = +el.value; fn(v); $(out).textContent = fmt(v); }; el.addEventListener('input', upd); upd(); };
  slider('#s3Iso', '#s3IsoV', v => uniforms.uIso.value = v, v => v.toFixed(2));
  slider('#s3Op', '#s3OpV', v => uniforms.uOpacity.value = v, v => v.toFixed(2));
  slider('#s3Speed', '#s3SpeedV', v => S.speed = v, v => v.toFixed(2) + '×');
  slider('#s3Steps', '#s3StepsV', v => uniforms.uSteps.value = v, v => v + ' pasos');

  /* ---------- operaciones ---------- */
  function runOp(op){
    if (S.anim) return;
    const a = SYM.analyze(S.g, S.mode, S.idx, S.acc, op);
    S.anim = { op, t: 0, a }; buildActive(op); describe(S.anim);
  }
  const OPDESC = { E: 'Identidad: no cambia nada.', C: 'Rotación propia de 360°/n alrededor del eje.', sigma: 'Reflexión en el plano: invierte la componente perpendicular. Se anima como un “paso a través” del plano: M(t) = I + (cos πt − 1) n̂n̂ᵀ.',
    i: 'Inversión r → −r. Se anima como S₂ = σₕ·C₂: rotación de 180° seguida de reflexión.', S: 'Rotación impropia Sₙ = σₕ·Cₙ: rota 360°/n y refleja en el plano perpendicular al eje.' };
  function matHTML(A){ return '<table class="s3mat">' + A.map(r => '<tr>' + r.map(x => `<td>${f3(x)}</td>`).join('') + '</tr>').join('') + '</table>'; }
  function termHTML(list){ return list.map((o, k) => { const c = o.c, s = c < 0 ? '−' : (k ? '+' : ''); const v = Math.abs(c); return `${s} ${Math.abs(v - 1) < 1e-6 ? '' : v.toFixed(3) + '·'}${SYM.BASIS[o.i]}`; }).join(' ') || '0'; }
  function describe(an){
    const G = SYM.GROUPS[S.g];
    const fname = S.mode === 'd' ? ORBN[S.idx] : 'SALC ' + G.salcs[S.idx].name.replace('_g', '<sub>g</sub>').replace('_u', '<sub>u</sub>');
    const irr = S.mode === 'd' ? G.dIrrep[S.idx] : G.salcs[S.idx].irrep;
    let h = `<div class="s3k">Función</div><div class="s3v">${fname} · representación ${G.irrepHTML[irr]} de ${G.name}</div>`;
    if (S.mode === 'salc') h += `<div class="s3note">Orbital del metal con la misma simetría: <b>${G.salcs[S.idx].mname}</b>. Solo funciones de la misma representación irreducible pueden combinarse para formar OM.</div>`;
    if (an) {
      const a = an.a, op = an.op;
      h += `<div class="s3k">Operación</div><div class="s3v">${op.label} · clase ${G.classes[op.cls]} · det = ${a.det > 0 ? '+1 (propia)' : '−1 (impropia)'}</div><div class="s3note">${OPDESC[op.type]}</div>`;
      h += `<div class="s3k">Antes → después</div><div class="s3v mono">${termHTML(a.before)}  →  ${termHTML(a.after)}</div>`;
      const ok = Math.abs(a.chiSet - a.table) < 1e-6;
      h += `<div class="s3chi"><div><span>χ(función)</span><b>${a.chi.toFixed(3).replace('-', '−')}</b><em>${a.kind}</em></div>
            <div><span>χ(${G.irrepHTML[a.irrep]}) calculado</span><b>${String(a.chiSet).replace('-', '−')}</b><em>tabla: ${String(a.table).replace('-', '−')} ${ok ? '✓' : '✗'}</em></div></div>`;
      $('#s3MatFinal').innerHTML = matHTML(a.R);
    } else { h += `<div class="s3note">Elige una operación: verás la matriz, la animación y el carácter calculado a partir de la función transformada.</div>`; $('#s3MatFinal').innerHTML = matHTML(SYM.I3()); }
    if (S.log.length) h += `<div class="s3k">Secuencia aplicada</div><div class="s3v mono">${S.log.join(' · ')}</div>`;
    $('#s3Info').innerHTML = h;
    // tabla de caracteres
    const ops = an ? an.op.cls : -1;
    let t = '<table class="s3ct"><thead><tr><th>' + G.name + '</th>' + G.classes.map((c, i) => `<th class="${i === ops ? 'on' : ''}">${c}</th>`).join('') + '</tr></thead><tbody>';
    Object.entries(G.irreps).forEach(([k, row]) => { t += `<tr class="${k === irr ? 'on' : ''}"><td>${G.irrepHTML[k]}</td>` + row.map((x, i) => `<td class="${i === ops ? 'on' : ''}">${String(x).replace('-', '−')}</td>`).join('') + '</tr>'; });
    $('#s3Table').innerHTML = t + '</tbody></table>';
  }

  /* ---------- bucle ---------- */
  function resize(){ const w = host.clientWidth, h = host.clientHeight; renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); }
  new ResizeObserver(() => { resize(); dirty = true; }).observe(host); resize();
  let visible = true, dirty = true, slowN = 0; new IntersectionObserver(es => { visible = es[0].isIntersecting; dirty = true; }).observe(host);
  controls.addEventListener('change', () => { dirty = true; });
  document.querySelector('.lab3d').addEventListener('input', () => { dirty = true; }, true); document.querySelector('.lab3d').addEventListener('click', () => { dirty = true; setTimeout(() => { dirty = true; }, 60); }, true);
  window.cristaliaLite = () => { visible = false; host.innerHTML = ''; sym3dLite(host); };
  const clock = new THREE.Clock(); let fpsT = 0, fpsN = 0, matTick = 0;
  const mtmp = new THREE.Matrix4();
  function frame(){
    requestAnimationFrame(frame);
    const dt = Math.min(clock.getDelta(), 0.05);
    if (!visible) return;
    const moved = controls.update();
    if (!S.anim && !moved && !dirty) return;                 // dibujar solo si algo cambió
    dirty = false; const t0f = performance.now();
    let M = S.acc;
    if (S.anim) {
      S.anim.t += dt * S.speed / 1.8;
      const t = S.anim.t >= 1 ? 1 : 0.5 - 0.5 * Math.cos(Math.PI * S.anim.t);          // suavizado
      M = SYM.mul(SYM.opMatrix(S.anim.op, t), S.acc);
      if (S.anim.t >= 1) { S.acc = SYM.clean(SYM.mul(SYM.finalMatrix(S.anim.op), S.acc)); S.log.push(S.anim.op.label); const done = S.anim; S.anim = null; M = S.acc; describe(done); }
    }
    const R3 = m3(M);
    mol.matrix.setFromMatrix3(R3); nodes.matrix.setFromMatrix3(R3);
    mol.matrixWorldNeedsUpdate = true; nodes.matrixWorldNeedsUpdate = true;   // necesario con matrixAutoUpdate = false
    uniforms.uMinv.value.copy(m3(SYM.inv(M)));
    camera.updateMatrixWorld();
    uniforms.uCam.value.copy(camera.position);
    uniforms.uViewProj.value.multiplyMatrices(camera.projectionMatrix, mtmp.copy(camera.matrixWorldInverse));
    renderer.render(scene, camera);
    // calidad adaptable: si el cuadro tarda mucho, bajar pasos y resolución
    const ft = performance.now() - t0f; slowN = ft > 28 ? slowN + 1 : Math.max(0, slowN - 1);
    if (slowN > 6 && uniforms.uSteps.value > 70) { uniforms.uSteps.value = Math.max(70, Math.round(uniforms.uSteps.value * 0.75)); const sl = $('#s3Steps'); if (sl) { sl.value = uniforms.uSteps.value; $('#s3StepsV').textContent = uniforms.uSteps.value + ' pasos (ajuste automático)'; } if (renderer.getPixelRatio() > 1) renderer.setPixelRatio(1); slowN = 0; dirty = true; }
    if (++matTick % 4 === 0) { $('#s3MatLive').innerHTML = matHTML(M) + `<div class="s3k" style="margin-top:4px">det = ${f3(SYM.det(M))}</div>`; }
    fpsN++; fpsT += dt; if (fpsT > 0.5) { $('#s3Fps').textContent = Math.round(fpsN / fpsT) + ' FPS'; fpsN = 0; fpsT = 0; }
  }
  resetAll(); frame();
}
