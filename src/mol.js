/* ===== Moléculas: Lewis, VSEPR, polaridad, OM; boranos, silanos, Hückel (puro) ===== */
const MOL = (function(){
"use strict";
const VAL={H:1,B:3,C:4,N:5,O:6,F:7,Si:4,P:5,S:6,Cl:7,Br:7,I:7,Xe:8,As:5};
const EN={H:2.20,B:2.04,C:2.55,N:3.04,O:3.44,F:3.98,Si:1.90,P:2.19,S:2.58,Cl:3.16,Br:2.96,I:2.66,Xe:2.60,As:2.18};
const nrm=v=>{const l=Math.hypot(...v)||1;return v.map(x=>x/l);};
const s3=1/Math.sqrt(3);
// direcciones VSEPR; los pares libres ocupan primero las posiciones indicadas
const DIRS={2:[[0,0,1],[0,0,-1]],3:[[0,1,0],[-.866,-.5,0],[.866,-.5,0]],4:[[s3,s3,s3],[-s3,-s3,s3],[-s3,s3,-s3],[s3,-s3,-s3]],
  5:[[1,0,0],[-.5,.866,0],[-.5,-.866,0],[0,0,1],[0,0,-1]],6:[[0,0,1],[0,0,-1],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]]};
const GEOMN={'2,0':'lineal','3,0':'trigonal plana','2,1':'angular','4,0':'tetraédrica','3,1':'piramidal trigonal','2,2':'angular','5,0':'bipirámide trigonal','4,1':'balancín','3,2':'forma de T','2,3':'lineal','6,0':'octaédrica','5,1':'pirámide cuadrada','4,2':'plano cuadrada','6,1':'octaédrica distorsionada (AX₆E)','1,0':'lineal (diatómica)','1,1':'lineal (diatómica)','1,2':'—','1,3':'—'};
const EDN={2:'lineal',3:'trigonal plana',4:'tetraédrica',5:'bipirámide trigonal',6:'octaédrica',7:'7 dominios (octaedro apuntado)'};
const HYB={2:'sp',3:'sp²',4:'sp³',5:'sp³d',6:'sp³d²',7:'sp³d³'};
// molécula de un átomo central: lig = [elemento, orden, pares libres del terminal]
function central(f,name,charge,c,lpC,lig,extra){
  const SN=lig.length+lpC, D=DIRS[SN]||DIRS[2]; const order=D.map((d,i)=>i);
  const lpSlots=order.slice(0,lpC), ligSlots=order.slice(lpC);
  const atoms=[{el:c,p:[0,0,0],lp:lpC}], bonds=[], lpDirs=lpSlots.map(i=>D[i]);
  lig.forEach(([el,o,lp],i)=>{ const L=el==='H'?0.75:1.0; const d=D[ligSlots[i]]; atoms.push({el,p:d.map(x=>x*L),lp}); bonds.push([0,i+1,o]); });
  return Object.assign({f,name,charge,atoms,bonds,lpDirs:{0:lpDirs}},extra||{});
}
const MOLS=[
 central('CH₄','metano',0,'C',0,[['H',1,0],['H',1,0],['H',1,0],['H',1,0]],{pg:'T<sub>d</sub>',org:true}),
 central('NH₃','amoníaco',0,'N',1,[['H',1,0],['H',1,0],['H',1,0]],{pg:'C<sub>3v</sub>'}),
 central('H₂O','agua',0,'O',2,[['H',1,0],['H',1,0]],{pg:'C<sub>2v</sub>'}),
 central('NH₄⁺','amonio',1,'N',0,[['H',1,0],['H',1,0],['H',1,0],['H',1,0]],{pg:'T<sub>d</sub>'}),
 central('CO₂','dióxido de carbono',0,'C',0,[['O',2,2],['O',2,2]],{pg:'D<sub>∞h</sub>',spec:'co2'}),
 central('NO₂⁺','ion nitronio',1,'N',0,[['O',2,2],['O',2,2]],{pg:'D<sub>∞h</sub>',note:'Isoelectrónico con CO₂ (16 e⁻ de valencia): lineal.'}),
 central('CO','monóxido de carbono',0,'C',1,[['O',3,1]],{pg:'C<sub>∞v</sub>',note:'Cargas formales C⁻/O⁺: el par libre del carbono es el que se une a los metales (σ-donador y π-aceptor).'}),
 central('HCN','cianuro de hidrógeno',0,'C',0,[['H',1,0],['N',3,1]],{pg:'C<sub>∞v</sub>',org:true}),
 central('O₃','ozono',0,'O',1,[['O',2,2],['O',1,3]],{pg:'C<sub>2v</sub>',note:'Dos formas resonantes equivalentes: los dos enlaces O–O son iguales (orden 1,5).'}),
 central('SO₂','dióxido de azufre',0,'S',1,[['O',2,2],['O',1,3]],{pg:'C<sub>2v</sub>',note:'Se muestra la estructura que cumple el octeto. La alternativa con dos S=O minimiza cargas formales, pero exige más de 8 e⁻ en S.'}),
 central('BF₃','trifluoruro de boro',0,'B',0,[['F',1,3],['F',1,3],['F',1,3]],{pg:'D<sub>3h</sub>',note:'B con solo 6 e⁻: ácido de Lewis; el orbital p vacío del B recibe algo de densidad π de los F.'}),
 central('NO₃⁻','nitrato',-1,'N',0,[['O',2,2],['O',1,3],['O',1,3]],{pg:'D<sub>3h</sub>',note:'Tres formas resonantes: los tres N–O son equivalentes.'}),
 central('CO₃²⁻','carbonato',-2,'C',0,[['O',2,2],['O',1,3],['O',1,3]],{pg:'D<sub>3h</sub>'}),
 central('SO₄²⁻','sulfato',-2,'S',0,[['O',1,3],['O',1,3],['O',1,3],['O',1,3]],{pg:'T<sub>d</sub>',note:'Estructura con octeto (S²⁺, O⁻). Los cálculos modernos favorecen esta descripción iónica frente a S=O con orbitales d.'}),
 central('XeO₆⁴⁻','ion perxenato',-4,'Xe',0,[['O',1,3],['O',1,3],['O',1,3],['O',1,3],['O',1,3],['O',1,3]],{pg:'O<sub>h</sub>',note:'Xe(VIII), octaédrico. La estructura con enlaces sencillos (Xe²⁺, O⁻) respeta el octeto de los O; las de Xe=O suponen expansión del octeto.'}),
 central('SF₆','hexafluoruro de azufre',0,'S',0,[['F',1,3],['F',1,3],['F',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'O<sub>h</sub>',note:'Molécula hipervalente: 12 e⁻ alrededor de S en la estructura de Lewis.'}),
 central('PCl₅','pentacloruro de fósforo',0,'P',0,[['Cl',1,3],['Cl',1,3],['Cl',1,3],['Cl',1,3],['Cl',1,3]],{pg:'D<sub>3h</sub>'}),
 central('SF₄','tetrafluoruro de azufre',0,'S',1,[['F',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>2v</sub>'}),
 central('ClF₃','trifluoruro de cloro',0,'Cl',2,[['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>2v</sub>'}),
 central('BrF₅','pentafluoruro de bromo',0,'Br',1,[['F',1,3],['F',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>4v</sub>'}),
 central('XeF₄','tetrafluoruro de xenón',0,'Xe',2,[['F',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'D<sub>4h</sub>'}),
 central('XeF₂','difluoruro de xenón',0,'Xe',3,[['F',1,3],['F',1,3]],{pg:'D<sub>∞h</sub>'}),
 central('I₃⁻','triyoduro',-1,'I',3,[['I',1,3],['I',1,3]],{pg:'D<sub>∞h</sub>'}),
];
// moléculas con varios centros (coordenadas explícitas)
(function(){
  // HClO₄: Cl tetraédrico; H sobre un O, fuera del eje
  const m=central('HClO₄','ácido perclórico',0,'Cl',0,[['O',1,3],['O',1,3],['O',1,3],['O',1,2]],{pg:'C<sub>s</sub>',note:'Estructura con octeto: Cl³⁺ y tres O⁻. La forma con tres Cl=O (cargas formales cero) requiere expansión del octeto; los cálculos actuales describen los enlaces como polares σ con interacción iónica.'});
  const o=m.atoms[4].p; m.atoms.push({el:'H',p:[o[0]+0.55,o[1]-0.55,o[2]+0.35],lp:0}); m.bonds.push([4,5,1]); MOLS.push(m);
  // B(OH)₃: plano, simetría C3h
  const at=[{el:'B',p:[0,0,0],lp:0}], bd=[];
  [90,210,330].forEach((a,i)=>{ const r=a*Math.PI/180; const O=[Math.cos(r),Math.sin(r),0]; at.push({el:'O',p:O,lp:2}); bd.push([0,at.length-1,1]);
    const rh=r+Math.PI/3; at.push({el:'H',p:[O[0]+0.7*Math.cos(rh),O[1]+0.7*Math.sin(rh),0],lp:0}); bd.push([at.length-2,at.length-1,1]); });
  MOLS.push({f:'B(OH)₃',name:'ácido bórico',charge:0,atoms:at,bonds:bd,pg:'C<sub>3h</sub>',note:'Ácido de Lewis (no de Brønsted): acepta OH⁻ del agua en el orbital p vacío del boro: B(OH)₃ + 2 H₂O ⇄ B(OH)₄⁻ + H₃O⁺.'});
  // C₂H₄
  MOLS.push({f:'C₂H₄',name:'eteno',charge:0,org:true,pg:'D<sub>2h</sub>',atoms:[{el:'C',p:[-.67,0,0],lp:0},{el:'C',p:[.67,0,0],lp:0},{el:'H',p:[-1.2,.9,0],lp:0},{el:'H',p:[-1.2,-.9,0],lp:0},{el:'H',p:[1.2,.9,0],lp:0},{el:'H',p:[1.2,-.9,0],lp:0}],bonds:[[0,1,2],[0,2,1],[0,3,1],[1,4,1],[1,5,1]]});
  // C₆H₆ (una estructura de Kekulé)
  const ca=[],cb=[]; for(let i=0;i<6;i++){ const a=Math.PI/2+i*Math.PI/3; ca.push({el:'C',p:[1.4*Math.cos(a),1.4*Math.sin(a),0],lp:0}); }
  for(let i=0;i<6;i++){ const a=Math.PI/2+i*Math.PI/3; ca.push({el:'H',p:[2.3*Math.cos(a),2.3*Math.sin(a),0],lp:0}); cb.push([i,(i+1)%6,i%2?1:2]); cb.push([i,i+6,1]); }
  MOLS.push({f:'C₆H₆',name:'benceno',charge:0,org:true,pg:'D<sub>6h</sub>',atoms:ca,bonds:cb,spec:'benzene',note:'Se dibuja una estructura de Kekulé; en realidad los seis C–C son iguales (1,39 Å): los 6 electrones π están deslocalizados.'});
})();

function analyzeMol(m){
  const n=m.atoms.length, deg=new Array(n).fill(0), bo=new Array(n).fill(0), pib=new Array(n).fill(0);
  m.bonds.forEach(([i,j,o])=>{ deg[i]++;deg[j]++; bo[i]+=o; bo[j]+=o; pib[i]+=o-1; pib[j]+=o-1; });
  const ve=m.atoms.reduce((s,a)=>s+VAL[a.el],0)-m.charge;
  const eLewis=2*m.bonds.reduce((s,b)=>s+b[2],0)+2*m.atoms.reduce((s,a)=>s+a.lp,0)+m.atoms.reduce((s,a)=>s+(a.rad||0),0);
  const fc=m.atoms.map((a,i)=>VAL[a.el]-2*a.lp-(a.rad||0)-bo[i]);
  const info=m.atoms.map((a,i)=>{ const SN=deg[i]+a.lp+(a.rad||0); return {SN,hyb:a.el==='H'?'—':(HYB[SN]||'—'),terminal:deg[i]<=1,geom:deg[i]>=2?GEOMN[deg[i]+','+(a.lp+(a.rad||0))]:null,edg:deg[i]>=2?EDN[SN]:null}; });
  // planaridad / linealidad
  const P=m.atoms.map(a=>a.p); const c0=P[0];
  let linear=true, planar=true;
  if(n>2){ let u=null; for(let i=1;i<n;i++){ const v=nrm(P[i].map((x,k)=>x-c0[k])); if(!u){u=v;continue;} const cr=Math.hypot(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]); if(cr>1e-3) linear=false; } }
  let normal=null;
  if(!linear){ outer: for(let i=1;i<n;i++)for(let j=i+1;j<n;j++){ const a=P[i].map((x,k)=>x-c0[k]), b=P[j].map((x,k)=>x-c0[k]); const cr=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]; if(Math.hypot(...cr)>1e-3){normal=nrm(cr);break outer;} }
    planar=P.every(p=>Math.abs(p.reduce((s,x,k)=>s+(x-c0[k])*normal[k],0))<1e-3); }
  // recuento de OM de valencia (base s,p de los átomos de los grupos principales)
  const ao=m.atoms.reduce((s,a)=>s+(a.el==='H'?1:4),0), somo=ve%2, occ=Math.ceil(ve/2), unocc=ao-occ;
  let piAO=0, piE=0, method;
  const heavy=m.atoms.map((a,i)=>({a,i})).filter(o=>o.a.el!=='H');
  const nPiBonds=m.bonds.reduce((s,b)=>s+(b[2]-1),0);
  if(linear&&n>=2){ method='simetría: eje molecular (dos sistemas π perpendiculares)'; piAO=2*heavy.length; piE=2*nPiBonds;
    heavy.forEach(({a,i})=>{ const axial=Math.max(0,2-deg[i]); const cap=2*(2-pib[i]); piE+=Math.min(cap,2*Math.max(0,a.lp-axial)); }); }
  else if(planar){ method='simetría: plano molecular (π = antisimétrico respecto al plano)'; piAO=heavy.length; piE=2*nPiBonds;
    heavy.forEach(({a,i})=>{ if(pib[i]===0&&a.lp>=1&&deg[i]+a.lp>=4) piE+=2; }); piE+=m.atoms.reduce((s,a)=>s+(a.rad||0),0); }
  else { method='enlace localizado (estructura de Lewis): sin plano ni eje, σ/π se asigna por enlace'; piAO=2*nPiBonds; piE=2*nPiBonds; }
  const piOcc=Math.ceil(piE/2), piUn=piAO-piOcc, sAO=ao-piAO, sOcc=occ-piOcc, sUn=sAO-sOcc;
  // momento dipolar por suma de dipolos de enlace (Pauling) + pares libres del centro (aprox.)
  let mu=[0,0,0];
  m.bonds.forEach(([i,j])=>{ const d=EN[m.atoms[j].el]-EN[m.atoms[i].el]; const u=nrm(P[j].map((x,k)=>x-P[i][k])); mu=mu.map((x,k)=>x+d*u[k]); });
  Object.entries(m.lpDirs||{}).forEach(([i,dirs])=>dirs.forEach(dv=>{ mu=mu.map((x,k)=>x+0.5*dv[k]); }));
  const muN=Math.hypot(...mu);
  return {ve,eLewis,fc,info,linear,planar,normal,ao,occ,unocc,somo,piAO,piE,piOcc,piUn,sAO,sOcc,sUn,method,mu,muN,polar:muN>0.05,deg};
}

/* ---------- espectro de masas: cúmulo isotópico del ion molecular ---------- */
const ISO={H:[[1.007825,99.9885],[2.014102,0.0115]],B:[[10.012937,19.9],[11.009305,80.1]],C:[[12,98.93],[13.003355,1.07]],N:[[14.003074,99.636],[15.000109,0.364]],
  O:[[15.994915,99.757],[16.999132,0.038],[17.99916,0.205]],F:[[18.998403,100]],Si:[[27.976927,92.223],[28.976495,4.685],[29.97377,3.092]],P:[[30.973762,100]],
  S:[[31.972071,94.99],[32.971459,0.75],[33.967867,4.25],[35.967081,0.01]],Cl:[[34.968853,75.76],[36.965903,24.24]],Br:[[78.918338,50.69],[80.916291,49.31]],I:[[126.904473,100]],As:[[74.921595,100]]};
function isotopePattern(m){
  if(m.atoms.some(a=>!ISO[a.el])) return null;
  let dist=new Map([[0,{p:1,mass:0}]]);
  m.atoms.forEach(a=>{ const nd=new Map(); dist.forEach((v,k)=>{ ISO[a.el].forEach(([ms,ab])=>{ const nk=k+Math.round(ms); const e=nd.get(nk)||{p:0,mass:0}; const pp=v.p*ab/100; e.mass=(e.mass*e.p+(v.mass+ms)*pp)/(e.p+pp); e.p+=pp; nd.set(nk,e); }); }); dist=new Map([...nd].filter(([k,v])=>v.p>1e-7)); });
  const arr=[...dist].map(([k,v])=>({mz:k,mass:v.mass,p:v.p})).sort((a,b)=>a.mz-b.mz); const mx=Math.max(...arr.map(x=>x.p));
  const mono=m.atoms.reduce((s,a)=>s+ISO[a.el][a.el==='B'?1:0][0],0);
  return {peaks:arr.map(x=>({mz:x.mz,mass:x.mass,rel:100*x.p/mx})).filter(x=>x.rel>0.05),mono};
}

/* ---------- Hückel ---------- */
function huckelRing(N,q){ const lev=[]; for(let k=0;k<N;k++) lev.push(2*Math.cos(2*Math.PI*k/N)); lev.sort((a,b)=>b-a); return fillHuckel(lev,N-q); } // energía en unidades de |β| (x: E = α + xβ)
/* ---------- constructor genérico (árbol VSEPR) y moléculas de los exámenes ---------- */
function rotTo(a,b){ // matriz que lleva el vector unitario a sobre b
  const v=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]], c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  if(c<-0.999999){ let p=Math.abs(a[0])<0.9?[1,0,0]:[0,1,0]; const k=nrm([a[1]*p[2]-a[2]*p[1],a[2]*p[0]-a[0]*p[2],a[0]*p[1]-a[1]*p[0]]); return axisRot(k,Math.PI); }
  const K=[[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]], I=[[1,0,0],[0,1,0],[0,0,1]], f=1/(1+c);
  return I.map((r,i)=>r.map((x,j)=>x+K[i][j]+f*[0,1,2].reduce((s,k)=>s+K[i][k]*K[k][j],0)));
}
function axisRot(k,t){ const c=Math.cos(t),s=Math.sin(t),C=1-c,[x,y,z]=k; return [[c+x*x*C,x*y*C-z*s,x*z*C+y*s],[y*x*C+z*s,c+y*y*C,y*z*C-x*s],[z*x*C-y*s,z*y*C+x*s,c+z*z*C]]; }
const mv=(M,v)=>M.map(r=>r[0]*v[0]+r[1]*v[1]+r[2]*v[2]);
// spec: [el, padre, orden, pares libres, giro(°), radical]
function tree(f,name,charge,spec,extra){
  const n=spec.length, atoms=spec.map(([el,,,lp,,rad])=>({el,lp,rad:rad||0,p:[0,0,0]})), bonds=[], nb=Array.from({length:n},()=>[]);
  spec.forEach(([,par,o],i)=>{ if(par>=0){ bonds.push([par,i,o]); nb[par].push(i); nb[i].push(par); } });
  const lpDirs={}, placed=new Array(n).fill(false); placed[0]=true;
  const place=i=>{ const a=atoms[i], SN=nb[i].length+a.lp+a.rad; const D=DIRS[SN]||DIRS[2]; const par=spec[i][1];
    let dirs;
    if(par<0) dirs=D.map(d=>d.slice());
    else { const u=nrm(a.p.map((x,k)=>x-atoms[par].p[k])); const toPar=u.map(x=>-x); const R=rotTo(nrm(D[0]),toPar); let ds=D.map(d=>mv(R,nrm(d)));
      // giro alrededor del enlace: anti respecto a otro vecino ya colocado del padre
      const other=nb[par].find(j=>j!==i&&placed[j]); if(ds.length>1){ let ref=other!=null?nrm(atoms[other].p.map((x,k)=>x-atoms[par].p[k])):null;
        if(ref){ const pr=v=>{ const d=v[0]*u[0]+v[1]*u[1]+v[2]*u[2]; return v.map((x,k)=>x-d*u[k]); }; const rp=pr(ref), dp=pr(ds[1]); if(Math.hypot(...rp)>1e-6&&Math.hypot(...dp)>1e-6){
          const r1=nrm(rp).map(x=>-x), d1=nrm(dp); let ang=Math.atan2([d1[1]*r1[2]-d1[2]*r1[1],d1[2]*r1[0]-d1[0]*r1[2],d1[0]*r1[1]-d1[1]*r1[0]].reduce((s,x,k)=>s+x*u[k],0),d1[0]*r1[0]+d1[1]*r1[1]+d1[2]*r1[2]);
          ang+=(spec[i][4]||0)*Math.PI/180; const Rz=axisRot(u,ang); ds=ds.map(d=>mv(Rz,d)); } } }
      dirs=ds; }
    const slots=dirs.slice(par<0?0:1); const lpS=slots.slice(0,a.lp+a.rad); if(lpS.length) lpDirs[i]=lpS.slice(0,a.lp);
    const kids=nb[i].filter(j=>j!==par); kids.forEach((j,q)=>{ const d=slots[a.lp+a.rad+q]||slots[q]; const L=(atoms[j].el==='H'||a.el==='H')?0.75:1.0; atoms[j].p=a.p.map((x,k)=>x+d[k]*L); placed[j]=true; });
    kids.forEach(place); };
  place(0);
  return Object.assign({f,name,charge,atoms,bonds,lpDirs},extra||{});
}
(function(){
  const s3=1/Math.sqrt(3), T=[[s3,s3,s3],[-s3,-s3,s3],[-s3,s3,-s3],[s3,-s3,-s3]];
  const add=m=>MOLS.push(m);
  add(central('SO₃²⁻','ion sulfito',-2,'S',1,[['O',1,3],['O',1,3],['O',1,3]],{pg:'C<sub>3v</sub>',note:'Estructura con octeto: S⁺ y tres O⁻ (carga total −2). Con enlaces S=O las cargas formales bajan, pero S supera el octeto.'}));
  add(central('SO₃','trióxido de azufre',0,'S',0,[['O',2,2],['O',1,3],['O',1,3]],{pg:'D<sub>3h</sub>',note:'Estructura con octeto (S²⁺); tres formas resonantes equivalentes: los tres S–O son iguales. La forma con tres S=O (cargas formales nulas) requiere 12 e⁻ alrededor de S.'}));
  add(central('XeO₄','tetróxido de xenón',0,'Xe',0,[['O',1,3],['O',1,3],['O',1,3],['O',1,3]],{pg:'T<sub>d</sub>',note:'Con octeto: Xe⁴⁺ y cuatro O⁻. La forma con cuatro Xe=O tiene cargas formales nulas.'}));
  add(central('XeO₃','trióxido de xenón',0,'Xe',1,[['O',1,3],['O',1,3],['O',1,3]],{pg:'C<sub>3v</sub>',note:'Pirámide trigonal (AX₃E), como NH₃. Explosivo.'}));
  add(central('IF₅','pentafluoruro de yodo',0,'I',1,[['F',1,3],['F',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>4v</sub>'}));
  add(central('BrF₃','trifluoruro de bromo',0,'Br',2,[['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>2v</sub>',note:'Forma de T (AX₃E₂): polar. Igual que ClF₃ (FCl₃ en el examen).'}));
  add(central('SOF₄','oxotetrafluoruro de azufre',0,'S',0,[['O',1,3],['F',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>2v</sub>',note:'Bipirámide trigonal (AX₅): el O ocupa una posición ecuatorial. Con octeto se escribe S⁺–O⁻.'}));
  add(central('SOCl₂','cloruro de tionilo',0,'S',1,[['O',1,3],['Cl',1,3],['Cl',1,3]],{pg:'C<sub>s</sub>',note:'Pirámide trigonal (AX₃E). Polar.'}));
  add(central('F₃PO','fluoruro de fosforilo',0,'P',0,[['O',1,3],['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>3v</sub>',note:'Tetraédrica. El enlace P–O es corto (polar, con carácter múltiple por interacción iónica), escrito P⁺–O⁻ con octeto o P=O con 10 e⁻.'}));
  add(central('ClO₄⁻','ion perclorato',-1,'Cl',0,[['O',1,3],['O',1,3],['O',1,3],['O',1,3]],{pg:'T<sub>d</sub>'}));
  add(central('CCl₄','tetracloruro de carbono',0,'C',0,[['Cl',1,3],['Cl',1,3],['Cl',1,3],['Cl',1,3]],{pg:'T<sub>d</sub>',note:'Enlaces C–Cl polares, pero la simetría tetraédrica anula el momento dipolar: apolar.'}));
  add(central('PH₃','fosfina',0,'P',1,[['H',1,0],['H',1,0],['H',1,0]],{pg:'C<sub>3v</sub>'}));
  add(central('CH₂Cl₂','diclorometano',0,'C',0,[['Cl',1,3],['H',1,0],['Cl',1,3],['H',1,0]],{pg:'C<sub>2v</sub>',org:true}));
  add(central('CHCl₃','cloroformo',0,'C',0,[['H',1,0],['Cl',1,3],['Cl',1,3],['Cl',1,3]],{pg:'C<sub>3v</sub>',org:true}));
  add(central('Cl₂CO','fosgeno',0,'C',0,[['O',2,2],['Cl',1,3],['Cl',1,3]],{pg:'C<sub>2v</sub>',org:true}));
  add(central('N₂O','óxido nitroso',0,'N',0,[['N',2,2],['O',2,2]],{pg:'C<sub>∞v</sub>',note:'N=N=O (dominante) y N≡N–O; lineal y ligeramente polar.'}));
  add(central('ClNO₂','cloruro de nitrilo',0,'N',0,[['Cl',1,3],['O',2,2],['O',1,3]],{pg:'C<sub>2v</sub>',note:'Plana; dos formas resonantes equivalentes para los N–O.'}));
  add((m=>{ m.atoms[0].rad=1; return m; })(central('ClO₂','dióxido de cloro',0,'Cl',1,[['O',2,2],['O',2,2]],{pg:'C<sub>2v</sub>',somoPi:true,note:'Radical con 19 e⁻ de valencia: un electrón desapareado (paramagnético). Se muestra la forma con cargas formales nulas (Cl con 2 enlaces dobles, 1 par libre y 1 electrón solo). El electrón desapareado ocupa un OM π* perpendicular al plano.'})));
  add(central('NF₃','trifluoruro de nitrógeno',0,'N',1,[['F',1,3],['F',1,3],['F',1,3]],{pg:'C<sub>3v</sub>',note:'AX₃E piramidal como el NH₃, pero casi apolar (μ = 0,23 D frente a 1,47 D): los dipolos de los enlaces N–F apuntan hacia los F y se oponen al del par libre del N. El par libre casi no es básico: el F se lo lleva.'}));
  (function(){ const m=central('NO₂','dióxido de nitrógeno',0,'N',1,[['O',2,2],['O',1,3]],{pg:'C<sub>2v</sub>',note:'Radical de 17 e⁻ de valencia (paramagnético, pardo rojizo). El electrón desapareado ocupa un orbital con mucho carácter del N, en el plano: por eso dos NO₂ se unen N–N formando N₂O₄ incoloro. El ángulo de 134° es mayor que 120°: un solo electrón repele menos que un par libre.'}); m.atoms[0].lp=0; m.atoms[0].rad=1; add(m); })();
  MOLS[MOLS.length-1].atoms[0].rad=1;
  add(central('C₂H₂*','(ver acetileno)',0,'C',0,[['C',3,0],['H',1,0]],{hide:true}));
  MOLS.pop();
  add(tree('C₂H₂','acetileno',0,[['C',-1,0,0],['C',0,3,0],['H',0,1,0],['H',1,1,0]],{pg:'D<sub>∞h</sub>',org:true}));
  add(tree('H₂CCO','cetena',0,[['C',-1,0,0],['C',0,2,0],['O',1,2,2],['H',0,1,0],['H',0,1,0]],{pg:'C<sub>2v</sub>',org:true,note:'CH₂ con C sp² y C central sp: plana. Los dos sistemas π (C=C y C=O) son perpendiculares entre sí.'}));
  add(tree('C₃H₈','propano',0,[['C',-1,0,0],['C',0,1,0],['C',0,1,0],['H',0,1,0],['H',0,1,0],['H',1,1,0],['H',1,1,0],['H',1,1,0],['H',2,1,0],['H',2,1,0],['H',2,1,0]],{pg:'C<sub>2v</sub>',org:true}));
  add(tree('H₂O₂','peróxido de hidrógeno',0,[['O',-1,0,2],['O',0,1,2,51.5],['H',0,1,0],['H',1,1,0]],{pg:'C<sub>2</sub>',note:'No plana: ángulo diedro H–O–O–H ≈ 111° en el gas. Polar.'}));
  add(tree('N₂H₄','hidrazina',0,[['N',-1,0,1],['N',0,1,1,-90],['H',0,1,0],['H',0,1,0],['H',1,1,0],['H',1,1,0]],{pg:'C<sub>2</sub>',note:'Conformación gauche: polar.'}));
  add(tree('N₂F₄','tetrafluoruro de dinitrógeno',0,[['N',-1,0,1],['N',0,1,1],['F',0,1,3],['F',0,1,3],['F',1,1,3],['F',1,1,3]],{pg:'C<sub>2h</sub>',note:'Existen las formas anti (C₂h, apolar) y gauche (polar) en equilibrio.'}));
  add(tree('N₂O₅','pentóxido de dinitrógeno',0,[['O',-1,0,2],['N',0,1,0],['N',0,1,0],['O',1,2,2],['O',1,1,3],['O',2,2,2],['O',2,1,3]],{pg:'C<sub>2</sub>',note:'O₂N–O–NO₂: cada N es trigonal plano (sp²) con carga formal +1; el sólido es iónico, NO₂⁺NO₃⁻.'}));
  add(tree('H₃PO₄','ácido fosfórico',0,[['P',-1,0,0],['O',0,1,3],['O',0,1,2],['O',0,1,2],['O',0,1,2],['H',2,1,0],['H',3,1,0],['H',4,1,0]],{pg:'C<sub>3v</sub>',note:'P tetraédrico; con octeto: P⁺–O⁻ (o P=O con 10 e⁻ y cargas nulas).'}));
  // P4
  { const P=T.map(v=>v.map(x=>x*1.3)), at=P.map(p=>({el:'P',p,lp:1})), b=[]; for(let i=0;i<4;i++) for(let j=i+1;j<4;j++) b.push([i,j,1]);
    add({f:'P₄',name:'fósforo blanco',charge:0,atoms:at,bonds:b,lpDirs:Object.fromEntries(T.map((v,i)=>[i,[v]])),pg:'T<sub>d</sub>',note:'Tetraedro: cada P con 3 enlaces y 1 par libre (sp³, piramidal). Ángulo P–P–P de 60°: enlaces muy tensionados, de ahí su alta reactividad.'}); }
  // As4O6 y As4O10
  const cage=(term)=>{ const A=T.map(v=>v.map(x=>x*1.55)), at=A.map(p=>({el:'As',p,lp:term?0:1})), b=[], lp={};
    for(let i=0;i<4;i++) for(let j=i+1;j<4;j++){ const mid=A[i].map((x,k)=>(x+A[j][k])/2), o=nrm(mid); at.push({el:'O',p:mid.map((x,k)=>x+o[k]*0.55),lp:2}); b.push([i,at.length-1,1],[j,at.length-1,1]); }
    if(term) T.forEach((v,i)=>{ at.push({el:'O',p:A[i].map((x,k)=>x+v[k]*1.0),lp:3}); b.push([i,at.length-1,1]); }); else T.forEach((v,i)=>lp[i]=[v]);
    return {atoms:at,bonds:b,lpDirs:lp}; };
  add(Object.assign({f:'As₄O₆',name:'óxido de arsénico(III)',charge:0,pg:'T<sub>d</sub>',note:'Jaula tipo adamantano: 4 As (piramidales, con un par libre) unidos por 6 O puente (angulares). Apolar.'},cage(false)));
  add(Object.assign({f:'As₄O₁₀',name:'óxido de arsénico(V)',charge:0,pg:'T<sub>d</sub>',note:'Como As₄O₆ pero cada As tiene además un O terminal: As tetraédrico (sp³). Con octeto: As⁺–O⁻ terminal. Apolar. Análogo a P₄O₁₀.'},cage(true)));
  // S3O9 cíclico
  { const at=[], b=[]; for(let i=0;i<6;i++){ const a=i*Math.PI/3; at.push({el:i%2?'O':'S',p:[1.45*Math.cos(a),1.45*Math.sin(a),i%2?0.25:-0.25],lp:i%2?2:0}); } for(let i=0;i<6;i++) b.push([i,(i+1)%6,1]);
    [0,2,4].forEach(i=>{ const p=at[i].p, o=nrm(p); [[o[0]*0.75,o[1]*0.75,0.66],[o[0]*0.75,o[1]*0.75,-0.66]].forEach(d=>{ at.push({el:'O',p:p.map((x,k)=>x+d[k]),lp:3}); b.push([i,at.length-1,1]); }); });
    add({f:'S₃O₉',name:'trímero cíclico del SO₃',charge:0,atoms:at,bonds:b,pg:'C<sub>3v</sub>',note:'Anillo de 6 miembros S–O–S–O–S–O: cada S tetraédrico con 2 O puente y 2 O terminales. Con octeto: S²⁺ y O⁻ terminales.'}); }
  // XeF6
  { const lpd=[s3,s3,s3], ax=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]].map(v=>nrm(v.map((x,k)=>x-0.22*lpd[k])));
    add({f:'XeF₆',name:'hexafluoruro de xenón',charge:0,atoms:[{el:'Xe',p:[0,0,0],lp:1}].concat(ax.map(v=>({el:'F',p:v.map(x=>x*1.1),lp:3}))),bonds:ax.map((_,i)=>[0,i+1,1]),lpDirs:{0:[lpd]},pg:'C<sub>3v</sub> (fluxional)',note:'AX₆E: 7 dominios. El par libre estereoactivo distorsiona el octaedro (octaedro apuntado, C₃ᵥ); la molécula es fluxional.'}); }
})();

function fillHuckel(lev,ne){ // lev ordenado de más estable (x mayor) a menos
  const occ=new Array(lev.length).fill(0); let left=ne; const groups=[];
  lev.forEach((x,i)=>{ const g=groups.find(g=>Math.abs(g.x-x)<1e-6); if(g) g.idx.push(i); else groups.push({x,idx:[i]}); });
  groups.forEach(g=>{ const cap=2*g.idx.length; const put=Math.min(cap,left); left-=put; let r=put; g.idx.forEach(i=>{ occ[i]=0; }); for(let pass=0;pass<2;pass++) g.idx.forEach(i=>{ if(r>0){ occ[i]++; r--; } }); });
  const E=lev.reduce((s,x,i)=>s+x*occ[i],0); let unp=occ.filter(o=>o===1).length;
  const homo=Math.max(...lev.map((x,i)=>occ[i]>0?-x:-Infinity)), lumo=Math.min(...lev.map((x,i)=>occ[i]<2?-x:Infinity));
  return {lev,occ,E,unp,gap:(lumo-homo)};
}
function huckelGraph(nAtoms,edges,ne){ const H=Array.from({length:nAtoms},()=>new Array(nAtoms).fill(0)); edges.forEach(([i,j])=>{H[i][j]=1;H[j][i]=1;}); const r=LF.jacobi(H); const lev=r.vals.slice().sort((a,b)=>b-a); return fillHuckel(lev,ne); }
function fusedRings(k){ // acenos lineales: k anillos
  const pts=[],edges=[]; const key=p=>p.map(x=>x.toFixed(3)).join(',');
  const idx=p=>{ let i=pts.findIndex(q=>key(q)===key(p)); if(i<0){pts.push(p);i=pts.length-1;} return i; };
  for(let r=0;r<k;r++){ const cx=r*Math.sqrt(3); const vs=[]; for(let i=0;i<6;i++){ const a=Math.PI/6+i*Math.PI/3; vs.push(idx([cx+Math.cos(a),Math.sin(a)])); }
    for(let i=0;i<6;i++){ const a=vs[i],b=vs[(i+1)%6]; if(!edges.some(e=>(e[0]===a&&e[1]===b)||(e[0]===b&&e[1]===a))) edges.push([a,b]); } }
  return {pts,edges};
}

/* ---------- boranos (styx de Lipscomb, reglas de Wade) ---------- */
function styxSolutions(m,n,c=0){ const out=[]; for(let s=0;s<=n-m;s++){ const x=n-m-s; const t=m-c-s; if(t<0) continue; const y=(3*m+n+c)/2-m-x-s-t; if(y<0||!Number.isInteger(y)) continue; out.push({s,t,y,x}); } return out; }
function wade(m,n,q){ const sep=(2*m+(n-m)+q)/2; const k=sep-m; return {sep,type:{1:'closo',2:'nido',3:'arachno',4:'hypho'}[k]||'—',parent:{1:m,2:m+1,3:m+2,4:m+3}[k]}; }

/* ---------- geometría experimental: longitudes (Å) y ángulos reales ---------- */
(function realGeometry(){
  const LEN={'C-H1':1.09,'H-N1':1.01,'H-O1':0.96,'H-P1':1.42,'B-O1':1.36,'B-F1':1.31,'C-O2':1.20,'C-O1':1.29,'C-O3':1.128,'C-N3':1.153,'C-Cl1':1.77,'C-C2':1.34,'C-C3':1.20,'C-C1':1.53,'N-O2':1.20,'N-O1':1.25,'N-N1':1.45,'N-N2':1.128,'F-N1':1.37,'Cl-N1':1.84,'O-O1':1.475,'O-O2':1.278,'Cl-P1':2.04,'O-P1':1.55,'F-P1':1.52,'F-S1':1.56,'O-S1':1.43,'O-S2':1.43,'Cl-S1':2.07,'Cl-F1':1.65,'Br-F1':1.75,'F-I1':1.86,'F-Xe1':1.95,'O-Xe1':1.76,'I-I1':2.92,'Cl-O1':1.44,'Cl-O2':1.47,'As-O1':1.79,'P-P1':2.21};
  // valores por molécula (fase gas salvo iones: valores típicos en sólidos)
  const MOLV={'CH₄':{},'NH₃':{ang:106.7},'H₂O':{ang:104.5,'H-O1':0.958},'NH₄⁺':{'H-N1':1.03},'CO₂':{'C-O2':1.162},'NO₂⁺':{'N-O2':1.154},'HCN':{'C-H1':1.065},
   'O₃':{ang:116.8,'O-O1':1.278,'O-O2':1.278},'SO₂':{ang:119.3},'NO₃⁻':{'N-O2':1.25,'N-O1':1.25},'CO₃²⁻':{'C-O2':1.29},'SO₄²⁻':{'O-S1':1.49,'O-S2':1.49},'XeO₆⁴⁻':{'O-Xe1':1.86},'SF₆':{'F-S1':1.561},
   'PCl₅':{ax:2.124,eq:2.02},'SF₄':{ang:101.6,ax:1.646,eq:1.545},'ClF₃':{ang:87.5,ax:1.698,eq:1.598},'BrF₅':{ang:84.8,ap:1.689,bas:1.774},'XeF₄':{'F-Xe1':1.95},'XeF₂':{'F-Xe1':1.977},
   'I₃⁻':{'I-I1':2.92},'SO₃²⁻':{ang:106,'O-S1':1.51,'O-S2':1.51},'SO₃':{'O-S1':1.42,'O-S2':1.42},'XeO₄':{'O-Xe1':1.736},'XeO₃':{ang:103},'IF₅':{ang:81.9,ap:1.84,bas:1.87},
   'BrF₃':{ang:86.2,ax:1.81,eq:1.72},'SOCl₂':{ang:96.1,'O-S1':1.44,'O-S2':1.44},'F₃PO':{'O-P1':1.44,'O-P2':1.44},'ClO₄⁻':{'Cl-O1':1.44,'Cl-O2':1.44},'PH₃':{ang:93.5},
   'Cl₂CO':{'C-O2':1.18,'C-Cl1':1.74},'ClO₂':{ang:117.4},'NF₃':{ang:102.4,'F-N1':1.365},'NO₂':{ang:134.3,'N-O2':1.197,'N-O1':1.197},'N₂O':{'N-O1':1.184,'N-O2':1.184},'ClNO₂':{'N-O2':1.20,'N-O1':1.20},'H₂CCO':{'C-C2':1.31,'C-O2':1.16},
   'XeF₆':{'F-Xe1':1.89},'C₆H₆':{'C-H1':1.08}};
  const key=(a,b,o)=>[a,b].sort().join('-')+o;
  const dist=(p,q)=>Math.hypot(p[0]-q[0],p[1]-q[1],p[2]-q[2]);
  const side=(m,i,j)=>{ // átomos del lado de j al cortar el enlace i-j (null si el enlace está en un anillo)
    const seen=new Set([j]), st=[j]; while(st.length){ const a=st.pop(); m.bonds.forEach(b=>{ let n=-1; if(b[0]===a) n=b[1]; else if(b[1]===a) n=b[0]; if(n<0||(a===j&&n===i)||seen.has(n)) return; if(n===i) return; seen.add(n); st.push(n); }); }
    return seen.has(i)?null:seen; };
  const inRing=(m,i,j)=>{ const seen=new Set([j]), st=[j]; while(st.length){ const a=st.pop(); for(const b of m.bonds){ let n=-1; if(b[0]===a) n=b[1]; else if(b[1]===a) n=b[0]; if(n<0||(a===j&&n===i)) continue; if(n===i) return true; if(!seen.has(n)){ seen.add(n); st.push(n); } } } return false; };
  MOLS.forEach(m=>{ const V=MOLV[m.f]||{}, A=m.atoms; m.expData=MOLV[m.f]||null;
    // 1) anillos y jaulas: escala uniforme según su enlace característico
    const rb=m.bonds.find(b=>inRing(m,b[0],b[1]));
    if(rb){ const want={'C₆H₆':1.397,'P₄':2.21,'As₄O₆':1.79,'As₄O₁₀':1.78,'S₃O₉':1.62}[m.f]; if(want){ const f=want/dist(A[rb[0]].p,A[rb[1]].p); A.forEach(a=>a.p=a.p.map(x=>x*f)); } }
    // 2) ángulos: los ligandos se alejan de los pares libres hasta el ángulo experimental
    const lp=(m.lpDirs&&m.lpDirs[0])||[]; const lig=m.bonds.filter(b=>b[0]===0).map(b=>b[1]);
    const cls={}; if(V.ax||V.ap){ const S0=lp.reduce((s,q)=>s.map((x,k)=>x+q[k]),[0,0,0]); lig.forEach(j=>{ const p=A[j].p, L=Math.hypot(...p), u=p.map(x=>x/L);
      if(V.ap){ const c=(u[0]*S0[0]+u[1]*S0[1]+u[2]*S0[2])/(Math.hypot(...S0)||1); cls[j]=c<-0.9?'ap':'bas'; }
      else cls[j]=(lig.some(k=>{ if(k===j) return false; const q=A[k].p, Lq=Math.hypot(...q); return (u[0]*q[0]+u[1]*q[1]+u[2]*q[2])/Lq<-0.95; })||lp.some(q=>u[0]*q[0]+u[1]*q[1]+u[2]*q[2]<-0.95))?'ax':'eq'; }); }
    if(m.f==='SF₄'){ const l=lp[0], nrm=v=>{ const L=Math.hypot(...v); return v.map(x=>x/L); }; const E=lig.filter(j=>cls[j]==='eq'), X=lig.filter(j=>cls[j]==='ax');
      const e1=nrm(A[E[0]].p), e2=nrm(A[E[1]].p), bis=nrm(e1.map((x,k)=>x+e2[k])), w=nrm(e1.map((x,k)=>x-e2[k])), h=101.6/2*Math.PI/180;
      [[E[0],1],[E[1],-1]].forEach(([j,sg])=>{ const L=Math.hypot(...A[j].p); A[j].p=bis.map((x,k)=>(Math.cos(h)*x+sg*Math.sin(h)*w[k])*L); });
      X.forEach(j=>{ const L=Math.hypot(...A[j].p), u=nrm(A[j].p), t=Math.tan((180-173.1)/2*Math.PI/180); A[j].p=nrm(u.map((x,k)=>x+t*bis[k])).map(x=>x*L); }); }
    else if(V.ang&&lp.length&&lig.length>=2){ const S=lp.reduce((s,d)=>s.map((x,k)=>x+d[k]),[0,0,0]); const D0=lig.map(i=>{ const p=A[i].p, L=Math.hypot(...p); return p.map(x=>x/L); });
      const nrm=v=>{ const L=Math.hypot(...v); return v.map(x=>x/L); }, angOf=(u,w)=>Math.acos(Math.max(-1,Math.min(1,u[0]*w[0]+u[1]*w[1]+u[2]*w[2])))*180/Math.PI;
      // par de referencia: para T, balancín y pirámide cuadrada se usa el ángulo indicado en MOLV.ref; si no, el menor entre ligandos
      const dirs=l=>D0.map(d=>{ const proj=d[0]*S[0]+d[1]*S[1]+d[2]*S[2]; const sl=Math.hypot(...S)||1; return nrm(d.map((x,k)=>x-l*S[k]/sl*(1+0*proj))); });
      const minAng=D=>{ let a=999; for(let x=0;x<D.length;x++)for(let y=x+1;y<D.length;y++) a=Math.min(a,angOf(D[x],D[y])); return a; };
      const SN=lig.length+lp.length; const target=V.ang;
      // para AX₄E (balancín) el ángulo de referencia es Feq–S–Feq (los dos ligandos más próximos al par libre)
      const refAng=D=>{ if(SN===5&&lig.length===4){ const idx=D0.map((d,i)=>[Math.abs(d[2]),i]).filter(x=>x[0]<0.5).map(x=>x[1]); return angOf(D[idx[0]],D[idx[1]]); } return minAng(D); };
      const f0=refAng(dirs(0)), dec=refAng(dirs(0.2))<f0; let lo=target<f0===dec?0:-0.42, hi=target<f0===dec?0.6:0;
      for(let it=0;it<60;it++){ const mid=(lo+hi)/2; const v=refAng(dirs(mid)); if((v>target)===dec) lo=mid; else hi=mid; }
      const Dn=dirs((lo+hi)/2); lig.forEach((i,k)=>{ const L=Math.hypot(...A[i].p); A[i].p=Dn[k].map(x=>x*L); });
      m.lpDirs[0]=lp.map(d=>nrm(d.map((x,k)=>x+0.15*d[k]))); }
    // 3) longitudes: cada enlace fuera de anillos al valor experimental (se mueve el fragmento terminal)
    m.bonds.forEach(b=>{ const [i,j,o]=b; if(inRing(m,i,j)) return; let want=V[key(A[i].el,A[j].el,o)]??LEN[key(A[i].el,A[j].el,o)]??(o>1?LEN[key(A[i].el,A[j].el,1)]:LEN[key(A[i].el,A[j].el,2)]);
      if(i===0&&cls[j]) want=V[cls[j]];
      if(!want) return; const sd=side(m,i,j); if(!sd) return; const L=dist(A[i].p,A[j].p); if(L<1e-6) return; const u=A[j].p.map((x,k)=>(x-A[i].p[k])/L), dl=want-L; sd.forEach(k=>{ A[k].p=A[k].p.map((x,q)=>x+u[q]*dl); }); });
    m.real=true; });
})();
return {MOLS,tree,analyzeMol,isotopePattern,huckelRing,huckelGraph,fusedRings,styxSolutions,wade,VAL,EN};
})();
if(typeof module!=='undefined') module.exports=MOL;
