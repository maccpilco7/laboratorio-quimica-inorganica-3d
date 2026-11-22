/* ===== Cristalografía, DRX y FRX (puro) ===== */
const CRY = (function(){
"use strict";
const D2R=Math.PI/180, NA=6.02214076e23, HC=12398.42; // eV·Å
// masas molares (g/mol, IUPAC) y Z
const EL={H:[1,1.008],B:[5,10.81],C:[6,12.011],N:[7,14.007],O:[8,15.999],F:[9,18.998],Na:[11,22.990],Mg:[12,24.305],Al:[13,26.982],Si:[14,28.085],P:[15,30.974],S:[16,32.06],Cl:[17,35.45],K:[19,39.098],Ca:[20,40.078],Ti:[22,47.867],V:[23,50.942],Cr:[24,51.996],Mn:[25,54.938],Fe:[26,55.845],Co:[27,58.933],Ni:[28,58.693],Cu:[29,63.546],Zn:[30,65.38],Sr:[38,87.62],Ag:[47,107.868],Cs:[55,132.905],Ba:[56,137.327],W:[74,183.84],Xe:[54,131.293],I:[53,126.904],Br:[35,79.904],Li:[3,6.94],Au:[79,196.967],Pb:[82,207.2],Zr:[40,91.224],'CO₂':[22,44.009],Po:[84,209],As:[33,74.922],Sb:[51,121.760],Sn:[50,118.710]};

/* ---------- sistemas cristalinos ---------- */
const SYSTEMS=[
 {k:'cub',n:'Cúbico',c:'a = b = c;  α = β = γ = 90°',cent:['P','I','F'],ex:'NaCl, Cu, CsCl, diamante',p:{a:4,b:4,c:4,al:90,be:90,ga:90},free:['a']},
 {k:'tet',n:'Tetragonal',c:'a = b ≠ c;  α = β = γ = 90°',cent:['P','I'],ex:'TiO₂ (rutilo), SnO₂',p:{a:4,b:4,c:6,al:90,be:90,ga:90},free:['a','c']},
 {k:'ort',n:'Ortorrómbico',c:'a ≠ b ≠ c;  α = β = γ = 90°',cent:['P','C','I','F'],ex:'marcasita FeS₂, azufre α, BaSO₄',p:{a:3.5,b:4.5,c:6,al:90,be:90,ga:90},free:['a','b','c']},
 {k:'hex',n:'Hexagonal',c:'a = b ≠ c;  α = β = 90°, γ = 120°',cent:['P'],ex:'Mg, Zr, ZnO, grafito, NiAs',p:{a:4,b:4,c:6.5,al:90,be:90,ga:120},free:['a','c']},
 {k:'tri',n:'Trigonal (romboédrico)',c:'a = b = c;  α = β = γ ≠ 90°',cent:['R'],ex:'antimonio, calcita, α-Al₂O₃ (corindón)',p:{a:4.5,b:4.5,c:4.5,al:70,be:70,ga:70},free:['a','al']},
 {k:'mon',n:'Monoclínico',c:'a ≠ b ≠ c;  α = γ = 90°, β ≠ 90°',cent:['P','C'],ex:'yeso, β-azufre',p:{a:4,b:5,c:6,al:90,be:110,ga:90},free:['a','b','c','be']},
 {k:'trc',n:'Triclínico',c:'a ≠ b ≠ c;  α ≠ β ≠ γ ≠ 90°',cent:['P'],ex:'CuSO₄·5H₂O, feldespatos',p:{a:4,b:5,c:6,al:80,be:105,ga:95},free:['a','b','c','al','be','ga']}
];
const CENT={P:[[0,0,0]],I:[[0,0,0],[.5,.5,.5]],F:[[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]],C:[[0,0,0],[.5,.5,0]],R:[[0,0,0]]};
const CENTN={P:'primitiva',I:'centrada en el cuerpo',F:'centrada en las caras',C:'centrada en las bases',R:'romboédrica (primitiva)'};
function cellVectors(p){ // matriz de vectores a, b, c (filas) en Å
  const al=p.al*D2R, be=p.be*D2R, ga=p.ga*D2R;
  const a=[p.a,0,0], b=[p.b*Math.cos(ga),p.b*Math.sin(ga),0];
  const cx=p.c*Math.cos(be), cy=p.c*(Math.cos(al)-Math.cos(be)*Math.cos(ga))/Math.sin(ga), cz=Math.sqrt(Math.max(0,p.c*p.c-cx*cx-cy*cy));
  return [a,b,[cx,cy,cz]];
}
function volume(p){ const ca=Math.cos(p.al*D2R),cb=Math.cos(p.be*D2R),cg=Math.cos(p.ga*D2R); return p.a*p.b*p.c*Math.sqrt(1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg); }
function frac2cart(f,M){ return [0,1,2].map(j=>f[0]*M[0][j]+f[1]*M[1][j]+f[2]*M[2][j]); }
// 1/d² = hᵀ G* h  (G* = inversa del tensor métrico)
function dspacing(p,h,k,l){
  const ca=Math.cos(p.al*D2R),cb=Math.cos(p.be*D2R),cg=Math.cos(p.ga*D2R);
  const G=[[p.a*p.a,p.a*p.b*cg,p.a*p.c*cb],[p.a*p.b*cg,p.b*p.b,p.b*p.c*ca],[p.a*p.c*cb,p.b*p.c*ca,p.c*p.c]];
  const det=G[0][0]*(G[1][1]*G[2][2]-G[1][2]*G[2][1])-G[0][1]*(G[1][0]*G[2][2]-G[1][2]*G[2][0])+G[0][2]*(G[1][0]*G[2][1]-G[1][1]*G[2][0]);
  const cof=(i,j)=>{const r=[0,1,2].filter(x=>x!==i),c=[0,1,2].filter(x=>x!==j);return ((i+j)%2?-1:1)*(G[r[0]][c[0]]*G[r[1]][c[1]]-G[r[0]][c[1]]*G[r[1]][c[0]]);};
  const Gi=[0,1,2].map(i=>[0,1,2].map(j=>cof(j,i)/det)); const v=[h,k,l];
  let s=0; for(let i=0;i<3;i++)for(let j=0;j<3;j++) s+=v[i]*Gi[i][j]*v[j];
  return 1/Math.sqrt(s);
}
// polígono del plano (hkl) dentro de la celda, en coordenadas fraccionarias
function planeAt(h,k,l,c){
  const V=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]], vals=V.map(v=>h*v[0]+k*v[1]+l*v[2]);
  const E=[[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]], pts=[];
  E.forEach(([i,j])=>{ const a=vals[i]-c,b=vals[j]-c; if(Math.abs(a)<1e-9) pts.push(V[i]); if(a*b<-1e-18){ const t=a/(a-b); pts.push(V[i].map((x,q)=>x+t*(V[j][q]-x))); } if(Math.abs(b)<1e-9) pts.push(V[j]); });
  const u0=[]; pts.forEach(p=>{ if(!u0.some(q=>Math.hypot(q[0]-p[0],q[1]-p[1],q[2]-p[2])<1e-6)) u0.push(p); });
  if(u0.length<3) return [];
  const cen=[0,1,2].map(q=>u0.reduce((s,p)=>s+p[q],0)/u0.length), n=[h,k,l];
  const u=[u0[0][0]-cen[0],u0[0][1]-cen[1],u0[0][2]-cen[2]], w=[n[1]*u[2]-n[2]*u[1],n[2]*u[0]-n[0]*u[2],n[0]*u[1]-n[1]*u[0]];
  const dot=(x,y)=>x[0]*y[0]+x[1]*y[1]+x[2]*y[2], ang=p=>{const pp=[p[0]-cen[0],p[1]-cen[1],p[2]-cen[2]]; return Math.atan2(dot(pp,w),dot(pp,u));};
  return u0.sort((p,q)=>ang(p)-ang(q));
}
function millerPolygon(h,k,l,shift=0){
  const V=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]];
  const vals=V.map(v=>h*v[0]+k*v[1]+l*v[2]); const mn=Math.min(...vals),mx=Math.max(...vals);
  let c=1+shift; if(!(c>mn+1e-9&&c<mx-1e-9)) c=(mn+mx)/2+shift;
  if(!(c>mn-1e-9&&c<mx+1e-9)) return [];
  const E=[[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]], pts=[];
  E.forEach(([i,j])=>{ const a=vals[i]-c,b=vals[j]-c; if(Math.abs(a)<1e-9) pts.push(V[i]); if(a*b<0){ const t=a/(a-b); pts.push(V[i].map((x,q)=>x+t*(V[j][q]-x))); } if(Math.abs(b)<1e-9) pts.push(V[j]); });
  const uniq=[]; pts.forEach(p=>{ if(!uniq.some(q=>Math.hypot(q[0]-p[0],q[1]-p[1],q[2]-p[2])<1e-6)) uniq.push(p); });
  if(uniq.length<3) return uniq;
  const cen=[0,1,2].map(q=>uniq.reduce((s,p)=>s+p[q],0)/uniq.length), n=[h,k,l];
  let u=[uniq[0][0]-cen[0],uniq[0][1]-cen[1],uniq[0][2]-cen[2]]; const w=[n[1]*u[2]-n[2]*u[1],n[2]*u[0]-n[0]*u[2],n[0]*u[1]-n[1]*u[0]];
  return uniq.sort((p,q)=>{ const pp=[p[0]-cen[0],p[1]-cen[1],p[2]-cen[2]],qq=[q[0]-cen[0],q[1]-cen[1],q[2]-cen[2]];
    const dot=(x,y)=>x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; return Math.atan2(dot(pp,w),dot(pp,u))-Math.atan2(dot(qq,w),dot(qq,u)); });
}

/* ---------- biblioteca de estructuras ---------- */
const fcc=[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]], add=(b,t)=>b.map(p=>[(p[0]+t[0])%1,(p[1]+t[1])%1,(p[2]+t[2])%1]);
const sp=(function(){ const A=[],B=[],O=[]; [[1/8,1/8,1/8],[7/8,3/8,3/8]].forEach(b=>fcc.forEach(t=>A.push(add([b],t)[0])));
  [[.5,.5,.5],[.25,.75,0],[.75,0,.25],[0,.25,.75]].forEach(b=>fcc.forEach(t=>B.push(add([b],t)[0])));
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)for(let k=0;k<4;k++){const o=(i%2)+(j%2)+(k%2); if(o===3||o===1) O.push([i/4,j/4,k/4]);} return {A,B,O}; })();
const STRUCT=[
 {k:'cu',n:'Cobre (fcc)',sys:'cub',p:{a:3.615},sites:[['Cu',fcc]],host:'fcc',note:'Metal cúbico centrado en las caras: empaquetamiento compacto ABC.'},
 {k:'al',n:'Aluminio (fcc)',sys:'cub',p:{a:4.0495},sites:[['Al',fcc]],host:'fcc'},
 {k:'fe',n:'Hierro α (bcc)',sys:'cub',p:{a:2.8665},sites:[['Fe',[[0,0,0],[.5,.5,.5]]]],host:'bcc',note:'Cúbico centrado en el cuerpo: NC = 8.'},
 {k:'w',n:'Wolframio (bcc)',sys:'cub',p:{a:3.165},sites:[['W',[[0,0,0],[.5,.5,.5]]]],host:'bcc'},
 {k:'nacl',n:'NaCl (sal gema)',sys:'cub',p:{a:5.640},sites:[['Na',fcc],['Cl',add(fcc,[.5,0,0])]],host:'fcc',note:'Cl⁻ en red fcc; Na⁺ ocupa todos los huecos octaédricos.'},
 {k:'kcl',n:'KCl (sal gema)',sys:'cub',p:{a:6.293},sites:[['K',fcc],['Cl',add(fcc,[.5,0,0])]],host:'fcc'},
 {k:'mgo',n:'MgO (periclasa)',sys:'cub',p:{a:4.212},sites:[['Mg',fcc],['O',add(fcc,[.5,0,0])]],host:'fcc'},
 {k:'nio',n:'NiO',sys:'cub',p:{a:4.177},sites:[['Ni',fcc],['O',add(fcc,[.5,0,0])]],host:'fcc'},
 {k:'cscl',n:'CsCl',sys:'cub',p:{a:4.123},sites:[['Cs',[[0,0,0]]],['Cl',[[.5,.5,.5]]]],note:'Dos redes cúbicas simples interpenetradas: NC = 8 : 8.'},
 {k:'zns',n:'ZnS (blenda)',sys:'cub',p:{a:5.409},sites:[['Zn',add(fcc,[.25,.25,.25])],['S',fcc]],host:'fcc',note:'S²⁻ fcc; Zn²⁺ en la mitad de los huecos tetraédricos.'},
 {k:'caf2',n:'CaF₂ (fluorita)',sys:'cub',p:{a:5.463},sites:[['Ca',fcc],['F',add(fcc,[.25,.25,.25]).concat(add(fcc,[.75,.75,.75]))]],host:'fcc',note:'Ca²⁺ fcc; F⁻ en todos los huecos tetraédricos: NC 8 : 4.'},
 {k:'si',n:'Silicio (diamante)',sys:'cub',p:{a:5.431},sites:[['Si',fcc.concat(add(fcc,[.25,.25,.25]))]],note:'Red fcc + la mitad de los huecos tetraédricos: cada átomo con NC = 4.'},
 {k:'dia',n:'Diamante (C)',sys:'cub',p:{a:3.567},sites:[['C',fcc.concat(add(fcc,[.25,.25,.25]))]]},
 {k:'srtio3',n:'SrTiO₃ (perovskita)',sys:'cub',p:{a:3.905},sites:[['Sr',[[0,0,0]]],['Ti',[[.5,.5,.5]]],['O',[[.5,.5,0],[.5,0,.5],[0,.5,.5]]]],note:'Ti⁴⁺ en octaedros TiO₆ unidos por vértices; Sr²⁺ con NC = 12.'},
 {k:'mgal2o4',n:'MgAl₂O₄ (espinela)',sys:'cub',p:{a:8.083},sites:[['Mg',sp.A],['Al',sp.B],['O',sp.O]],note:'Espinela normal (u ideal).'},
 {k:'fe3o4',n:'Fe₃O₄ (magnetita)',sys:'cub',p:{a:8.396},sites:[['Fe',sp.A.concat(sp.B)],['O',sp.O]],note:'Espinela inversa (posiciones de u ideal).'},
 {k:'tio2',n:'TiO₂ (rutilo)',sys:'tet',p:{a:4.594,c:2.959},sites:[['Ti',[[0,0,0],[.5,.5,.5]]],['O',[[.305,.305,0],[.695,.695,0],[.805,.195,.5],[.195,.805,.5]]]],note:'Octaedros TiO₆ que comparten aristas en cadenas a lo largo de c.'},
 {k:'zno',n:'ZnO (wurtzita)',sys:'hex',p:{a:3.250,c:5.207},sites:[['Zn',[[1/3,2/3,0],[2/3,1/3,.5]]],['O',[[1/3,2/3,.382],[2/3,1/3,.882]]]],note:'O²⁻ hcp; Zn²⁺ en la mitad de los huecos tetraédricos.'},
 {k:'po',n:'Polonio α (cúbica simple)',sys:'cub',p:{a:3.359},sites:[['Po',[[0,0,0]]]],note:'Único elemento con red cúbica simple: NC = 6 y solo 52 % del volumen ocupado.'},
 {k:'au',n:'Oro (fcc)',sys:'cub',p:{a:4.078},sites:[['Au',fcc]],host:'fcc',note:'Examen: con ρ = 19,3 g/cm³ y Z = 4 se obtiene a ≈ 4,08 Å.'},
 {k:'lif',n:'LiF (sal gema)',sys:'cub',p:{a:4.03},sites:[['Li',add(fcc,[.5,0,0])],['F',fcc]],host:'fcc',note:'F⁻ en red fcc y Li⁺ en todos los huecos octaédricos. Sus planos (200), d = 2,014 Å, se usan como cristal analizador LiF(200) en FRX.'},
 {k:'licl',n:'LiCl (sal gema)',sys:'cub',p:{a:5.14},sites:[['Li',add(fcc,[.5,0,0])],['Cl',fcc]],host:'fcc',note:'Arista del examen: 5,14 Å.'},
 {k:'kbr',n:'KBr (sal gema)',sys:'cub',p:{a:6.578},sites:[['K',add(fcc,[.5,0,0])],['Br',fcc]],host:'fcc',note:'Arista del examen: 6,578 Å.'},
 {k:'agcl',n:'AgCl (sal gema)',sys:'cub',p:{a:5.549},sites:[['Ag',add(fcc,[.5,0,0])],['Cl',fcc]],host:'fcc',note:'Con ρ = 5,56 g/cm³: Z = ρ·N_A·a³/M ≈ 4 → celda centrada en las caras.'},
 {k:'pbs',n:'PbS (galena)',sys:'cub',p:{a:5.936},sites:[['Pb',add(fcc,[.5,0,0])],['S',fcc]],host:'fcc',note:'Estructura tipo NaCl.'},
 {k:'co2',n:'CO₂ sólido (hielo seco)',sys:'cub',p:{a:5.540},sites:[['CO₂',fcc]],host:'fcc',note:'Moléculas de CO₂ en los nudos de una red fcc (cada esfera representa una molécula). Arista del examen: 5,540 Å.'},
 {k:'mno2',n:'MnO₂ (pirolusita)',sys:'tet',p:{a:4.39,c:2.86},sites:[['Mn',[[0,0,0],[.5,.5,.5]]],['O',[[.305,.305,0],[.695,.695,0],[.805,.195,.5],[.195,.805,.5]]]],note:'Tipo rutilo (examen: a = 4,39 Å, c = 2,86 Å). Z = 2.'},
 {k:'zr',n:'Circonio (hcp)',sys:'hex',p:{a:3.232,c:5.147},sites:[['Zr',[[1/3,2/3,.25],[2/3,1/3,.75]]]],host:'hcp',note:'Examen: a = 323,2 pm, c = 514,7 pm.'},
 {k:'mg',n:'Magnesio (hcp)',sys:'hex',p:{a:3.209,c:5.211},sites:[['Mg',[[1/3,2/3,.25],[2/3,1/3,.75]]]],host:'hcp',note:'Empaquetamiento compacto hexagonal ABAB: NC = 12.'}
 ,{k:'grafito',n:'Grafito (hexagonal)',sys:'hex',p:{a:2.464,c:6.711},sites:[['C',[[0,0,.25],[0,0,.75],[1/3,2/3,.25],[2/3,1/3,.75]]]],note:'Capas de hexágonos (C–C 1,42 Å, sp²) apiladas ABAB a 3,35 Å; entre capas solo fuerzas de van der Waals.'}
 ,{k:'nias',n:'NiAs (arseniuro de níquel)',sys:'hex',p:{a:3.619,c:5.034},sites:[['Ni',[[0,0,0],[0,0,.5]]],['As',[[1/3,2/3,.25],[2/3,1/3,.75]]]],note:'As en hcp; Ni en todos los huecos octaédricos. Ni: NC 6 (octaedro); As: NC 6 (prisma trigonal). Tipo de muchos sulfuros y arseniuros metálicos.'}
 ,{k:'sb',n:'Antimonio (trigonal, ejes hexagonales)',sys:'hex',p:{a:4.307,c:11.273},sites:[['Sb',(()=>{ const u=0.2336, B=[[0,0,u],[0,0,1-u]], T=[[0,0,0],[2/3,1/3,1/3],[1/3,2/3,2/3]]; const o=[]; T.forEach(t=>B.forEach(b=>o.push([(b[0]+t[0])%1,(b[1]+t[1])%1,(b[2]+t[2])%1]))); return o; })()]],note:'Sistema trigonal (red romboédrica R descrita con ejes hexagonales, 6 átomos por celda). Capas plegadas: cada Sb tiene 3 vecinos cercanos y 3 más lejanos (3 + 3). Mineral de origen: estibina, Sb₂S₃, abundante en Bolivia.'}
 ,{k:'casiterita',n:'Casiterita SnO₂ (rutilo)',sys:'tet',p:{a:4.737,c:3.186},sites:[['Sn',[[0,0,0],[.5,.5,.5]]],['O',[[.307,.307,0],[.693,.693,0],[.807,.193,.5],[.193,.807,.5]]]],note:'Mena principal del estaño (Huanuni, Colquiri, Llallagua). Estructura de rutilo: Sn⁴⁺ en octaedros SnO₆ (NC 6) y O con NC 3. Por FRX se ve el Sn; por DRX se distingue de otras fases con Sn.'}
 ,{k:'marcasita',n:'Marcasita FeS₂ (ortorrómbica)',sys:'ort',p:{a:4.443,b:5.423,c:3.387},sites:[['Fe',[[0,0,0],[.5,.5,.5]]],['S',[[.2002,.378,0],[.7998,.622,0],[.2998,.878,.5],[.7002,.122,.5]]]],note:'Ortorrómbica (Pnnm). Fe en octaedros FeS₆; el azufre forma pares S₂²⁻ (S–S ≈ 2,2 Å), como en la pirita cúbica. Mineral común en yacimientos de sulfuros.'}

];
function fullParams(s){ const S=SYSTEMS.find(x=>x.k===s.sys); const p=Object.assign({},S.p,s.p); if(s.sys==='cub'){p.b=p.c=p.a;} if(s.sys==='tet'||s.sys==='hex'){p.b=p.a;} return p; }
// átomos de la celda (posiciones únicas en [0,1))
function cellAtoms(s){ const out=[]; s.sites.forEach(([el,list])=>list.forEach(f=>{ const g=f.map(x=>((x%1)+1)%1); if(!out.some(a=>a.el===el&&a.f.every((x,i)=>Math.abs(x-g[i])<1e-6))) out.push({el,f:g}); })); return out; }
// átomos de dibujo (incluye imágenes en los bordes) con su peso en la celda
function drawAtoms(s){ const out=[]; cellAtoms(s).forEach(a=>{ const opts=a.f.map(x=>Math.abs(x)<1e-6?[0,1]:[x]);
  opts[0].forEach(x=>opts[1].forEach(y=>opts[2].forEach(z=>{ const nb=[x,y,z].filter(v=>Math.abs(v)<1e-6||Math.abs(v-1)<1e-6).length; out.push({el:a.el,f:[x,y,z],w:[1,.5,.25,.125][nb],pos:['interior','cara','arista','vértice'][nb]}); })));
  }); return out; }
function gcd(a,b){ return b?gcd(b,a%b):a; }
function analyze(s){
  const p=fullParams(s), V=volume(p), M=cellVectors(p), at=cellAtoms(s), dr=drawAtoms(s);
  const count={}; at.forEach(a=>count[a.el]=(count[a.el]||0)+1);
  const byPos={}; dr.forEach(a=>{ const k=a.el+'|'+a.pos; byPos[k]=byPos[k]||{el:a.el,pos:a.pos,n:0,w:a.w}; byPos[k].n++; });
  const els=Object.keys(count); let g=0; els.forEach(e=>g=gcd(g,count[e])); const Zf=g; const formula=els.map(e=>({el:e,n:count[e]/g}));
  const Mform=formula.reduce((s,f)=>s+f.n*EL[f.el][1],0), rho=Zf*Mform/(NA*V*1e-24);
  // números de coordinación: por cada especie, vecinos más cercanos de cada especie contraria (o de sí misma si es un elemento)
  const ANION=['O','F','Cl','S','Br','I','N'], cn={};
  const near=(a,b)=>{ let best=Infinity; for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)for(let k=-1;k<=1;k++){ const c=frac2cart([b.f[0]+i-a.f[0],b.f[1]+j-a.f[1],b.f[2]+k-a.f[2]],M); const d=Math.hypot(...c); if(d>1e-6) best=Math.min(best,d); } return best; };
  const countAt=(a,list,d0)=>{ let n=0; list.forEach(b=>{ for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)for(let k=-1;k<=1;k++){ const c=frac2cart([b.f[0]+i-a.f[0],b.f[1]+j-a.f[1],b.f[2]+k-a.f[2]],M); const d=Math.hypot(...c); if(d>1e-6&&d<d0*1.03) n++; } }); return n; };
  const isAnion=e=>ANION.includes(e)&&els.length>1, partnersOf=e=>els.length===1?[e]:els.filter(f=>f!==e&&(isAnion(e)?!ANION.includes(f):ANION.includes(f)||!els.some(g=>ANION.includes(g))));
  const dist=(a,b)=>{ let best=Infinity; for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)for(let k=-1;k<=1;k++){ const c=frac2cart([b.f[0]+i-a.f[0],b.f[1]+j-a.f[1],b.f[2]+k-a.f[2]],M); const d=Math.hypot(...c); if(d>1e-6&&d<best) best=d; } return best; };
  const shell=(a,list,dl)=>{ let n=0; list.forEach(b=>{ for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)for(let k=-1;k<=1;k++){ const c=frac2cart([b.f[0]+i-a.f[0],b.f[1]+j-a.f[1],b.f[2]+k-a.f[2]],M); const d=Math.hypot(...c); if(d>1e-6&&d<dl(b)) n++; } }); return n; };
  // primera capa de cada átomo (para elementos en varios sitios, como Fe en Fe3O4)
  const d0=at.map(a=>Math.min(...at.filter(b=>partnersOf(a.el).includes(b.el)).map(b=>dist(a,b))));
  at.forEach((a,i)=>a._d0=d0[i]);
  els.forEach(e=>{ const pe=partnersOf(e), mine=at.filter(x=>x.el===e);
    const siteOf=a=>{ const parts=pe.map(f=>{ const list=at.filter(x=>x.el===f); const n=isAnion(e)?shell(a,list,b=>b._d0*1.03):shell(a,list,()=>a._d0*1.03); return {el:f,n,d:Math.min(...list.map(b=>dist(a,b)))}; }); return {n:parts.reduce((s,x)=>s+x.n,0),parts,d:a._d0}; };
    const sites=[]; mine.forEach(a=>{ const S=siteOf(a); const k=sites.find(t=>t.n===S.n&&Math.abs(t.d-S.d)<0.02); if(k) k.m++; else sites.push(Object.assign(S,{m:1})); });
    const dmax=Math.max(...sites.map(t=>t.d)), dmin=Math.min(...sites.map(t=>t.d));
    const parts=pe.map(f=>({el:f,n:sites.length===1?sites[0].parts.find(p=>p.el===f).n:sites.map(t=>t.parts.find(p=>p.el===f).n).join(' o '),d:dmax}));
    cn[e]={n:sites.length===1?sites[0].n:sites.map(t=>t.n).join(' y '),parts,d:dmin,dmax,sites,multi:sites.length>1}; });
  return {p,V,M,count,byPos:Object.values(byPos),formula,Z:Zf,Mform,rho,cn,at,dr};
}
// huecos de una red fcc (posiciones fraccionarias, sin duplicar)
function holes(type){
  if(type==='fcc') return {oct:[[.5,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5]],tet:[[.25,.25,.25],[.75,.25,.25],[.25,.75,.25],[.25,.25,.75],[.75,.75,.25],[.75,.25,.75],[.25,.75,.75],[.75,.75,.75]]};
  if(type==='bcc') return {oct:[[.5,.5,0],[.5,0,.5],[0,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5]],tet:[[.5,.25,0],[.5,.75,0],[.25,.5,0],[.75,.5,0],[.5,0,.25],[.5,0,.75],[.25,0,.5],[.75,0,.5],[0,.5,.25],[0,.5,.75],[0,.25,.5],[0,.75,.5]]};
  if(type==='hcp') return {oct:[[0,0,0],[0,0,.5]],tet:[[1/3,2/3,.625],[1/3,2/3,.875],[2/3,1/3,.125],[2/3,1/3,.375]]};
  return null;
}

/* ---------- difractograma de rayos X (polvo) ---------- */
// f(s) aproximado: Z·exp(−4,5 s²) con s = senθ/λ (aproximación didáctica, sin términos de Cromer–Mann)
function xrd(s,lambda,tmax=120){
  const p=fullParams(s), at=cellAtoms(s), peaks={};
  const lim=7;
  for(let h=-lim;h<=lim;h++)for(let k=-lim;k<=lim;k++)for(let l=-lim;l<=lim;l++){
    if(!h&&!k&&!l) continue; const d=dspacing(p,h,k,l); const sth=lambda/(2*d); if(sth>=1) continue; const th=Math.asin(sth); const tt=2*th/D2R; if(tt>tmax) continue;
    const ss=sth/lambda; let re=0,im=0; at.forEach(a=>{ const f=EL[a.el][0]*Math.exp(-4.5*ss*ss); const ph=2*Math.PI*(h*a.f[0]+k*a.f[1]+l*a.f[2]); re+=f*Math.cos(ph); im+=f*Math.sin(ph); });
    const F2=re*re+im*im; if(F2<1e-3) continue;
    const LP=(1+Math.cos(2*th)**2)/(Math.sin(th)**2*Math.cos(th));
    const key=tt.toFixed(3); if(!peaks[key]) peaks[key]={tt,d,I:0,hkl:[h,k,l],mult:0};
    peaks[key].I+=F2*LP; peaks[key].mult++;
    const cur=peaks[key].hkl; const score=v=>(v[0]>=0&&v[1]>=0&&v[2]>=0?0:10)+v.filter(x=>x<0).length; if(score([h,k,l])<score(cur)||(score([h,k,l])===score(cur)&&(h>cur[0]||(h===cur[0]&&k>cur[1])||(h===cur[0]&&k===cur[1]&&l>cur[2])))) peaks[key].hkl=[h,k,l];
  }
  const list=Object.values(peaks).sort((a,b)=>a.tt-b.tt); const mx=Math.max(...list.map(p=>p.I));
  list.forEach(p=>p.rel=100*p.I/mx); return list.filter(p=>p.rel>=0.5);
}
const XLINES={ // energías en eV — X-Ray Data Booklet, tabla 1-2 (K…Zn) y Kaye & Laby (Ag)
  K:[3313.8,3311.1,3589.6],Ca:[3691.68,3688.09,4012.7],Sc:[4090.6,4086.1,4460.5],Ti:[4510.84,4504.86,4931.81],V:[4952.20,4944.64,5427.29],Cr:[5414.72,5405.509,5946.71],
  Mn:[5898.75,5887.65,6490.45],Fe:[6403.84,6390.84,7057.98],Co:[6930.32,6915.30,7649.43],Ni:[7478.15,7460.89,8264.66],Cu:[8047.78,8027.83,8905.29],Zn:[8638.86,8615.78,9572.0],Ag:[22163,21991,24943],
  Na:[1040.98,1040.98,1071.1],Mg:[1253.60,1253.60,1302.2],Al:[1486.70,1486.27,1557.45],Si:[1739.98,1739.38,1835.94],P:[2013.7,2012.7,2139.1],S:[2307.84,2306.64,2464.04],Cl:[2622.39,2620.78,2815.6],
  Br:[11924.2,11877.6,13291.4],Sr:[14165,14097.9,15835.7],I:[28612.0,28317.2,32294.7],Cs:[30972.8,30625.1,34986.9],Ba:[32193.6,31817.1,36378.2],W:[59318.24,57981.7,67244.3]};
const TUBES={Cu:12398.42/8047.78,Co:12398.42/6930.32,Fe:12398.42/6403.84,Cr:12398.42/5414.72};
// líneas K de un elemento: tabla verificada o estimación de Moseley (marcada)
function kLines(el){ const Z=EL[el][0]; if(XLINES[el]){ const [a1,a2,b]=XLINES[el]; return {Z,ka1:a1,ka2:a2,kb:b,est:false}; }
  if(Z<11) return null; const ka=10.2*(Z-1)*(Z-1), kb=12.09*(Z-1.8)*(Z-1.8); return {Z,ka1:ka,ka2:ka-1.5e-4*ka,kb,est:true}; }
const CRYSTALS=[{k:'LiF200',n:'LiF(200)',d2:4.028},{k:'LiF220',n:'LiF(220)',d2:2.848},{k:'PG002',n:'Grafito pirolítico (002)',d2:6.715},{k:'PET',n:'PET(002)',d2:8.742},{k:'TAP',n:'TAP(001)',d2:25.7}];
function twoTheta(E,d2){ const lam=HC/E; const s=lam/d2; return s<1?2*Math.asin(s)/D2R:null; }
const omegaK=Z=>{ const z4=Math.pow(Z,4); return z4/(1.12e6+z4); }; // rendimiento de fluorescencia K (aproximación tipo Burhop)
return {planeAt,SYSTEMS,CENT,CENTN,cellVectors,volume,frac2cart,dspacing,millerPolygon,STRUCT,fullParams,cellAtoms,drawAtoms,analyze,holes,xrd,XLINES,TUBES,kLines,CRYSTALS,twoTheta,omegaK,EL,HC};
})();
if(typeof module!=='undefined') module.exports=CRY;
