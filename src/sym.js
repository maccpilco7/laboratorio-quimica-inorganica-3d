/* ===== Motor de simetría molecular (puro, sin dependencias) =====
   Convención química: z vertical. Matrices 3×3 como arrays de filas.
   Representación activa: (O_R f)(p) = f(R⁻¹ p).                              */
const SYM = (function(){
"use strict";
const I3=()=>[[1,0,0],[0,1,0],[0,0,1]];
const mul=(A,B)=>A.map((r,i)=>[0,1,2].map(j=>r[0]*B[0][j]+r[1]*B[1][j]+r[2]*B[2][j]));
const mv=(A,v)=>A.map(r=>r[0]*v[0]+r[1]*v[1]+r[2]*v[2]);
const tr=A=>[0,1,2].map(i=>[0,1,2].map(j=>A[j][i]));
const det=A=>A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
function inv(A){ const d=det(A); const c=(i,j)=>{const r=[0,1,2].filter(x=>x!==i), s=[0,1,2].filter(x=>x!==j); return A[r[0]][s[0]]*A[r[1]][s[1]]-A[r[0]][s[1]]*A[r[1]][s[0]];};
  return [0,1,2].map(i=>[0,1,2].map(j=>((i+j)%2?-1:1)*c(j,i)/d)); }
const nrm=v=>{const l=Math.hypot(...v);return v.map(x=>x/l);};
function rot(n,th){ n=nrm(n); const [x,y,z]=n, c=Math.cos(th), s=Math.sin(th), C=1-c;
  return [[c+x*x*C, x*y*C-z*s, x*z*C+y*s],[y*x*C+z*s, c+y*y*C, y*z*C-x*s],[z*x*C-y*s, z*y*C+x*s, c+z*z*C]]; }
function squash(n,c){ n=nrm(n); return [0,1,2].map(i=>[0,1,2].map(j=>(i===j?1:0)+(c-1)*n[i]*n[j])); } // c=−1 → reflexión
const clean=A=>A.map(r=>r.map(x=>Math.abs(x)<1e-12?0:Math.abs(x-Math.round(x))<1e-12?Math.round(x):x));

/* ---- operación interpolada M(t), t∈[0,1] ---- */
function opMatrix(op,t){
  const s=Math.max(0,Math.min(1,t));
  const safe=c=>Math.abs(c)<0.035?(c<0?-0.035:0.035):c;
  if(op.type==='E') return I3();
  if(op.type==='C') return rot(op.n, s*2*Math.PI/op.k);
  if(op.type==='sigma') return squash(op.n, s>=1?-1:safe(Math.cos(Math.PI*s)));
  if(op.type==='S'||op.type==='i'){ // S_n = σh·C_n ; i = S₂
    const k=op.type==='i'?2:op.k, n=op.type==='i'?[0,0,1]:op.n, th=2*Math.PI/k;
    if(s<0.5) return rot(n, 2*s*th);
    const u=2*s-1; return mul(squash(n, u>=1?-1:safe(Math.cos(Math.PI*u))), rot(n,th));
  }
}
const finalMatrix=op=>clean(opMatrix(op,1));

/* ---- base de funciones: 5 d + 3 p + s + 8 ligandos ---- */
const r3=Math.sqrt(3)/2;
const AM=[ [[-.5,0,0],[0,-.5,0],[0,0,1]], [[r3,0,0],[0,-r3,0],[0,0,0]], [[0,r3,0],[r3,0,0],[0,0,0]], [[0,0,r3],[0,0,0],[r3,0,0]], [[0,0,0],[0,0,r3],[0,r3,0]] ];
const frob=(A,B)=>{let s=0;for(let i=0;i<3;i++)for(let j=0;j<3;j++)s+=A[i][j]*B[i][j];return s;};
const BASIS=['dz²','dx²−y²','dxy','dxz','dyz','px','py','pz','s','L1','L2','L3','L4','L5','L6','L7','L8'];
for(let k=1;k<=8;k++) for(const a of 'xyz') BASIS.push('L'+k+'·p'+a);   // orbitales p de los ligandos (enlace π)
const NB=41;
// T(R): f' = T f  con (O_R f)(p) = f(Rᵀp), R ortogonal
function Tmat(R, ligs){
  const T=Array.from({length:NB},()=>new Array(NB).fill(0));
  for(let i=0;i<5;i++){ const RA=mul(mul(R,AM[i]),tr(R)); for(let j=0;j<5;j++) T[j][i]=frob(AM[j],RA)/1.5; }
  for(let i=0;i<3;i++) for(let j=0;j<3;j++) T[5+j][5+i]=R[j][i];
  T[8][8]=1;
  ligs.forEach((L,k)=>{ const RL=mv(R,L); const j=ligs.findIndex(M=>Math.hypot(M[0]-RL[0],M[1]-RL[1],M[2]-RL[2])<1e-6); if(j>=0){ T[9+j][9+k]=1;
    for(let a=0;a<3;a++) for(let b=0;b<3;b++) T[17+3*j+b][17+3*k+a]=R[b][a]; } });   // vector p en L_k → R·p en L_j
  return T;
}
const matvec=(T,f)=>T.map(r=>r.reduce((s,x,i)=>s+x*f[i],0));
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
function decompose(f,tol=1e-6){ const n=Math.sqrt(dot(f,f)); return f.map((c,i)=>({c:c/n,i})).filter(o=>Math.abs(o.c)>tol); }

/* ---- grupos ---- */
const s3=1/Math.sqrt(3), s2=1/Math.SQRT2;
const D=2.6; // distancia M–L en unidades de la escena
const OHL=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
const TDL=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(v=>v.map(x=>x*s3));
const SQL=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]];
const pad=a=>a.concat(new Array(8-a.length).fill(0));
const GROUPS={
  Oh:{name:'O<sub>h</sub>',geom:'Octaedro ML₆',ligs:OHL,
    classes:['E','8C₃','6C₂′','6C₄','3C₂ (=C₄²)','i','6S₄','8S₆','3σₕ','6σd'],
    ops:[{id:'E',label:'E',type:'E',cls:0},{id:'C4',label:'C₄ (z)',type:'C',n:[0,0,1],k:4,cls:3},{id:'C2',label:'C₂ (z)',type:'C',n:[0,0,1],k:2,cls:4},
         {id:'C3',label:'C₃ [111]',type:'C',n:[1,1,1],k:3,cls:1},{id:'C2p',label:'C₂′ [110]',type:'C',n:[1,1,0],k:2,cls:2},{id:'i',label:'i',type:'i',cls:5},
         {id:'S4',label:'S₄ (z)',type:'S',n:[0,0,1],k:4,cls:6},{id:'S6',label:'S₆ [111]',type:'S',n:[1,1,1],k:6,cls:7},
         {id:'sh',label:'σₕ (xy)',type:'sigma',n:[0,0,1],cls:8},{id:'sd',label:'σd (⊥[1−10])',type:'sigma',n:[1,-1,0],cls:9}],
    irreps:{A1g:[1,1,1,1,1,1,1,1,1,1],Eg:[2,-1,0,0,2,2,0,-1,2,0],T2g:[3,0,1,-1,-1,3,-1,0,-1,1],T1u:[3,0,-1,1,-1,-3,-1,0,1,1]},
    irrepHTML:{A1g:'a<sub>1g</sub>',Eg:'e<sub>g</sub>',T2g:'t<sub>2g</sub>',T1u:'t<sub>1u</sub>'},
    dIrrep:['Eg','Eg','T2g','T2g','T2g'], dSets:{Eg:[0,1],T2g:[2,3,4]},
    salcs:[{name:'a₁g',irrep:'A1g',set:'a1g',c:[1,1,1,1,1,1].map(x=>x/Math.sqrt(6)),metal:8,mname:'4s'},
           {name:'e_g (z²)',irrep:'Eg',set:'eg',c:[-1,-1,-1,-1,2,2].map(x=>x/Math.sqrt(12)),metal:0,mname:'d<sub>z²</sub>'},
           {name:'e_g (x²−y²)',irrep:'Eg',set:'eg',c:[1,1,-1,-1,0,0].map(x=>x/2),metal:1,mname:'d<sub>x²−y²</sub>'},
           {name:'t₁u (x)',irrep:'T1u',set:'t1u',c:[1,-1,0,0,0,0].map(x=>x*s2),metal:5,mname:'4p<sub>x</sub>'},
           {name:'t₁u (y)',irrep:'T1u',set:'t1u',c:[0,0,1,-1,0,0].map(x=>x*s2),metal:6,mname:'4p<sub>y</sub>'},
           {name:'t₁u (z)',irrep:'T1u',set:'t1u',c:[0,0,0,0,1,-1].map(x=>x*s2),metal:7,mname:'4p<sub>z</sub>'},
           {name:'π t₂g (xy)',irrep:'T2g',set:'t2gpi',pi:[[0,1,0],[0,-1,0],[1,0,0],[-1,0,0],[0,0,0],[0,0,0]].map(v=>v.map(x=>x/2)),metal:2,mname:'d<sub>xy</sub>'},
           {name:'π t₂g (xz)',irrep:'T2g',set:'t2gpi',pi:[[0,0,1],[0,0,-1],[0,0,0],[0,0,0],[1,0,0],[-1,0,0]].map(v=>v.map(x=>x/2)),metal:3,mname:'d<sub>xz</sub>'},
           {name:'π t₂g (yz)',irrep:'T2g',set:'t2gpi',pi:[[0,0,0],[0,0,0],[0,0,1],[0,0,-1],[0,1,0],[0,-1,0]].map(v=>v.map(x=>x/2)),metal:4,mname:'d<sub>yz</sub>'}],
    elements:()=>{const e=[];[[1,0,0],[0,1,0],[0,0,1]].forEach(n=>e.push({t:'C',n,k:4}));[[1,1,1],[1,1,-1],[1,-1,1],[-1,1,1]].forEach(n=>e.push({t:'C',n,k:3}));
      [[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]].forEach(n=>e.push({t:'C',n,k:2,p:1}));[[1,0,0],[0,1,0],[0,0,1]].forEach(n=>e.push({t:'sigma',n,kind:'h'}));
      [[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]].forEach(n=>e.push({t:'sigma',n,kind:'d'}));e.push({t:'i'});return e;},
    edges:[[0,2],[0,3],[0,4],[0,5],[1,2],[1,3],[1,4],[1,5],[2,4],[2,5],[3,4],[3,5]]},
  Td:{name:'T<sub>d</sub>',geom:'Tetraedro ML₄',ligs:TDL,
    classes:['E','8C₃','3C₂','6S₄','6σd'],
    ops:[{id:'E',label:'E',type:'E',cls:0},{id:'C3',label:'C₃ [111]',type:'C',n:[1,1,1],k:3,cls:1},{id:'C2',label:'C₂ (z)',type:'C',n:[0,0,1],k:2,cls:2},
         {id:'S4',label:'S₄ (z)',type:'S',n:[0,0,1],k:4,cls:3},{id:'sd',label:'σd (⊥[1−10])',type:'sigma',n:[1,-1,0],cls:4}],
    irreps:{A1:[1,1,1,1,1],E:[2,-1,2,0,0],T2:[3,0,-1,-1,1]},
    irrepHTML:{A1:'a<sub>1</sub>',E:'e',T2:'t<sub>2</sub>'},
    dIrrep:['E','E','T2','T2','T2'], dSets:{E:[0,1],T2:[2,3,4]},
    salcs:[{name:'a₁',irrep:'A1',set:'a1',c:[1,1,1,1].map(x=>x/2),metal:8,mname:'4s'},
           {name:'t₂ (x)',irrep:'T2',set:'t2',c:[1,1,-1,-1].map(x=>x/2),metal:5,mname:'4p<sub>x</sub> (o d<sub>yz</sub>)'},
           {name:'t₂ (y)',irrep:'T2',set:'t2',c:[1,-1,1,-1].map(x=>x/2),metal:6,mname:'4p<sub>y</sub> (o d<sub>xz</sub>)'},
           {name:'t₂ (z)',irrep:'T2',set:'t2',c:[1,-1,-1,1].map(x=>x/2),metal:7,mname:'4p<sub>z</sub> (o d<sub>xy</sub>)'}],
    elements:()=>{const e=[];[[1,1,1],[1,1,-1],[1,-1,1],[-1,1,1]].forEach(n=>e.push({t:'C',n,k:3}));[[1,0,0],[0,1,0],[0,0,1]].forEach(n=>e.push({t:'C',n,k:2,s4:1}));
      [[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]].forEach(n=>e.push({t:'sigma',n,kind:'d'}));return e;},
    edges:[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]},
  D4h:{name:'D<sub>4h</sub>',geom:'Plano cuadrado ML₄',ligs:SQL,
    classes:['E','2C₄','C₂','2C₂′','2C₂″','i','2S₄','σₕ','2σᵥ','2σd'],
    ops:[{id:'E',label:'E',type:'E',cls:0},{id:'C4',label:'C₄ (z)',type:'C',n:[0,0,1],k:4,cls:1},{id:'C2',label:'C₂ (z)',type:'C',n:[0,0,1],k:2,cls:2},
         {id:'C2p',label:'C₂′ (x)',type:'C',n:[1,0,0],k:2,cls:3},{id:'C2pp',label:'C₂″ [110]',type:'C',n:[1,1,0],k:2,cls:4},{id:'i',label:'i',type:'i',cls:5},
         {id:'S4',label:'S₄ (z)',type:'S',n:[0,0,1],k:4,cls:6},{id:'sh',label:'σₕ (xy)',type:'sigma',n:[0,0,1],cls:7},
         {id:'sv',label:'σᵥ (xz)',type:'sigma',n:[0,1,0],cls:8},{id:'sd',label:'σd (⊥[1−10])',type:'sigma',n:[1,-1,0],cls:9}],
    irreps:{A1g:[1,1,1,1,1,1,1,1,1,1],B1g:[1,-1,1,1,-1,1,-1,1,1,-1],B2g:[1,-1,1,-1,1,1,-1,1,-1,1],Eg:[2,0,-2,0,0,2,0,-2,0,0],Eu:[2,0,-2,0,0,-2,0,2,0,0]},
    irrepHTML:{A1g:'a<sub>1g</sub>',B1g:'b<sub>1g</sub>',B2g:'b<sub>2g</sub>',Eg:'e<sub>g</sub>',Eu:'e<sub>u</sub>'},
    dIrrep:['A1g','B1g','B2g','Eg','Eg'], dSets:{A1g:[0],B1g:[1],B2g:[2],Eg:[3,4]},
    salcs:[{name:'a₁g',irrep:'A1g',set:'a1g',c:[1,1,1,1].map(x=>x/2),metal:0,mname:'d<sub>z²</sub> (y 4s)'},
           {name:'b₁g',irrep:'B1g',set:'b1g',c:[1,1,-1,-1].map(x=>x/2),metal:1,mname:'d<sub>x²−y²</sub>'},
           {name:'e_u (x)',irrep:'Eu',set:'eu',c:[1,-1,0,0].map(x=>x*s2),metal:5,mname:'4p<sub>x</sub>'},
           {name:'e_u (y)',irrep:'Eu',set:'eu',c:[0,0,1,-1].map(x=>x*s2),metal:6,mname:'4p<sub>y</sub>'}],
    elements:()=>{const e=[{t:'C',n:[0,0,1],k:4}];[[1,0,0],[0,1,0]].forEach(n=>e.push({t:'C',n,k:2,p:1}));[[1,1,0],[1,-1,0]].forEach(n=>e.push({t:'C',n,k:2,p:2}));
      e.push({t:'sigma',n:[0,0,1],kind:'h'});[[1,0,0],[0,1,0]].forEach(n=>e.push({t:'sigma',n,kind:'v'}));[[1,1,0],[1,-1,0]].forEach(n=>e.push({t:'sigma',n,kind:'d'}));e.push({t:'i'});return e;},
    edges:[[0,2],[2,1],[1,3],[3,0]]}
};
// vector de función en la base combinada
function fnVector(g,mode,idx,overlay){
  const f=new Array(NB).fill(0);
  if(mode==='d'){ f[idx]=1; return f; }
  const sl=GROUPS[g].salcs[idx]; if(sl.c) pad(sl.c).forEach((c,k)=>f[9+k]=c);
  if(sl.pi) sl.pi.forEach((v,k)=>{ for(let a=0;a<3;a++) f[17+3*k+a]=v[a]; });
  return f;
}
function renderVector(g,mode,idx,overlay){
  const f=fnVector(g,mode,idx); if(mode==='salc'&&overlay){ const sl=GROUPS[g].salcs[idx]; f[sl.metal]+=0.85; } return f;
}
function setVectors(g,mode,idx){
  const G=GROUPS[g];
  if(mode==='d'){ const ir=G.dIrrep[idx]; return {irrep:ir, vecs:G.dSets[ir].map(i=>fnVector(g,'d',i))}; }
  const sl=G.salcs[idx]; return {irrep:sl.irrep, vecs:G.salcs.map((s,i)=>({s,i})).filter(o=>o.s.set===sl.set).map(o=>fnVector(g,'salc',o.i))};
}
// análisis de una operación R sobre la función actual (orientación acumulada Acc)
function analyze(g,mode,idx,Acc,op){
  const G=GROUPS[g], ligs=G.ligs, R=finalMatrix(op);
  const TA=Tmat(Acc,ligs), TR=Tmat(R,ligs);
  const f=matvec(TA,fnVector(g,mode,idx)), Rf=matvec(TR,f);
  const chi=dot(f,Rf)/dot(f,f);
  const set=setVectors(g,mode,idx); let chiSet=0; set.vecs.forEach(v=>{const a=matvec(TA,v); chiSet+=dot(a,matvec(TR,a))/dot(a,a);});
  const table=G.irreps[set.irrep][op.cls];
  let kind; const rnd=x=>Math.round(x*1e6)/1e6;
  if(Math.abs(chi-1)<1e-6) kind='se mantiene (+1)'; else if(Math.abs(chi+1)<1e-6) kind='invierte su fase (−1)'; else if(Math.abs(chi)<1e-6 && decompose(Rf).length===1) kind='pasa a otra función del conjunto (χ = 0)'; else kind='se transforma en una combinación lineal';
  return {R, det:det(R), f, Rf, chi:rnd(chi), chiSet:rnd(chiSet), table, irrep:set.irrep, kind, before:decompose(f), after:decompose(Rf)};
}
return {I3,mul,mv,tr,det,inv,rot,opMatrix,finalMatrix,Tmat,matvec,dot,decompose,GROUPS,BASIS,NB,D,renderVector,fnVector,setVectors,analyze,clean};
})();
if(typeof module!=='undefined') module.exports=SYM;

/* ===== Laboratorio de simetría sin conexión (Canvas 2D, sin Three.js) ===== */
function sym3dLite(host){
  const G=SYM.GROUPS; let g='Oh', mode='d', idx=0, Acc=SYM.I3(), anim=null, iso=0.3, yaw=0.6, pitch=0.45, drag=null, pts=[];
  host.innerHTML=''; host.style.background='#0b1018';
  const cv=document.createElement('canvas'); cv.style.cssText='width:100%;height:100%;display:block;touch-action:none;cursor:grab'; host.appendChild(cv);
  const panel=document.querySelector('.lab3d .panel'); const note=document.createElement('div'); note.className='s3msg'; note.style.cssText='position:absolute;left:10px;top:8px;font-size:12px;color:#9fb0c8;background:none;padding:0';
  note.textContent='Modo sin conexión: vista simplificada (nube de puntos) sin WebGL ni Three.js.'; host.appendChild(note);
  const AM=[[[-.5,0,0],[0,-.5,0],[0,0,1]],[[.866,0,0],[0,-.866,0],[0,0,0]],[[0,.866,0],[.866,0,0],[0,0,0]],[[0,0,.866],[0,0,0],[.866,0,0]],[[0,0,0],[0,0,.866],[0,.866,0]]];
  const psi=(f,x,y,z,L)=>{ const r=Math.hypot(x,y,z), e=Math.exp(-r*1.1); let v=0; const q=[x,y,z];
    for(let k=0;k<5;k++) if(f[k]){ const A=AM[k]; let s=0; for(let i=0;i<3;i++)for(let j=0;j<3;j++) s+=q[i]*A[i][j]*q[j]; v+=f[k]*s*e; }
    for(let a=0;a<3;a++) if(f[5+a]) v+=f[5+a]*q[a]*e*1.4; if(f[8]) v+=f[8]*e*0.8;
    L.forEach((l,k)=>{ const dx=x-l[0]*2.2, dy=y-l[1]*2.2, dz=z-l[2]*2.2, g2=Math.exp(-(dx*dx+dy*dy+dz*dz)/0.35); if(f[9+k]) v+=f[9+k]*g2*1.2; for(let a=0;a<3;a++) if(f[17+3*k+a]) v+=f[17+3*k+a]*[dx,dy,dz][a]*g2*2; }); return v; };
  function sample(){ const f=SYM.renderVector(g,mode,idx,false), L=G[g].ligs; pts=[]; const N=34, R=3.2; let mx=0; const vals=[];
    for(let i=0;i<N;i++)for(let j=0;j<N;j++)for(let k=0;k<N;k++){ const x=-R+2*R*i/(N-1), y=-R+2*R*j/(N-1), z=-R+2*R*k/(N-1); const v=psi(f,x,y,z,L); vals.push([x,y,z,v]); mx=Math.max(mx,Math.abs(v)); }
    const th=iso*mx; vals.forEach(([x,y,z,v])=>{ const a=Math.abs(v); if(a>th&&a<th*1.35) pts.push([x+(Math.random()-.5)*0.12,y+(Math.random()-.5)*0.12,z+(Math.random()-.5)*0.12,v>0?1:-1]); }); }
  const mul=SYM.mul;
  function draw(){ const dpr=window.devicePixelRatio||1, W=cv.clientWidth||600, H=cv.clientHeight||500; cv.width=W*dpr; cv.height=H*dpr; const ctx=cv.getContext('2d'); ctx.setTransform(dpr,0,0,dpr,0,0); ctx.clearRect(0,0,W,H);
    const M=anim?mul(SYM.opMatrix(anim.op,anim.t),Acc):Acc; const ca=Math.cos(yaw),sa=Math.sin(yaw),cb=Math.cos(pitch),sb=Math.sin(pitch), S=Math.min(W,H)/7.5;
    const P=p=>{ const q=SYM.mv(M,p); const x1=q[0]*ca-q[1]*sa, y1=q[0]*sa+q[1]*ca, y2=y1*cb-q[2]*sb, z2=y1*sb+q[2]*cb; return [W/2+x1*S,H/2-z2*S,y2]; };
    const Pf=p=>{ const x1=p[0]*ca-p[1]*sa, y1=p[0]*sa+p[1]*ca, y2=y1*cb-p[2]*sb, z2=y1*sb+p[2]*cb; return [W/2+x1*S,H/2-z2*S,y2]; };
    const L=G[g].ligs.map(l=>l.map(x=>x*2.2));
    if(document.querySelector('#s3Poly')&&document.querySelector('#s3Poly').checked) (G[g].edges||[]).forEach(([a,b])=>{ const p=P(L[a]), q=P(L[b]); ctx.strokeStyle='rgba(160,180,210,.45)'; ctx.lineWidth=1.2; ctx.beginPath(); ctx.moveTo(p[0],p[1]); ctx.lineTo(q[0],q[1]); ctx.stroke(); });
    const items=pts.map(p=>{ const q=P(p); return [q[0],q[1],q[2],p[3]]; }).sort((a,b)=>b[2]-a[2]);
    items.forEach(([x,y,z,s])=>{ const al=Math.max(.25,Math.min(.95,.6-z/8)); ctx.fillStyle=s>0?`rgba(79,134,255,${al})`:`rgba(255,138,61,${al})`; ctx.fillRect(x-1.6,y-1.6,3.2,3.2); });
    L.forEach((l,k)=>{ const q=P(l); ctx.fillStyle=k===0?'#f2b84b':'#9fb0c8'; ctx.beginPath(); ctx.arc(q[0],q[1],7,0,7); ctx.fill(); ctx.fillStyle='#0b1018'; ctx.font='600 9px IBM Plex Sans, sans-serif'; ctx.textAlign='center'; ctx.fillText('L'+(k+1),q[0],q[1]+3); ctx.textAlign='left'; });
    const c=P([0,0,0]); ctx.fillStyle='#d8dee9'; ctx.beginPath(); ctx.arc(c[0],c[1],6,0,7); ctx.fill();
    if(anim&&anim.op.n){ const a=Pf(anim.op.n.map(x=>x*3.4)), b=Pf(anim.op.n.map(x=>-x*3.4)); ctx.strokeStyle='#f2b84b'; ctx.setLineDash([6,4]); ctx.lineWidth=2; ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke(); ctx.setLineDash([]); } }
  function buildUI(){ panel.innerHTML=`<div class="ctl"><label>Grupo puntual</label><div class="seg" id="l3G">${Object.keys(G).map(k=>`<button data-k="${k}" aria-pressed="${k===g}">${G[k].name||k}</button>`).join('')}</div></div>
    <div class="ctl"><label>Función</label><div class="seg" id="l3F">${SYM.BASIS.slice(0,5).map((n,i)=>`<button data-m="d" data-i="${i}" aria-pressed="${mode==='d'&&i===idx}">${n}</button>`).join('')}${G[g].salcs.map((s,i)=>`<button data-m="salc" data-i="${i}" aria-pressed="${mode==='salc'&&i===idx}">SALC ${s.name}</button>`).join('')}</div></div>
    <div class="ctl"><label>Operaciones (se encadenan)</label><div class="ops" id="l3O">${G[g].ops.map((o,i)=>`<button data-o="${i}">${o.label}</button>`).join('')}<button data-o="r">Reiniciar</button></div></div>
    <div class="ctl"><label for="l3I">Isovalor: <b id="l3IV">${iso}</b></label><input type="range" id="l3I" min="0.1" max="0.8" step="0.05" value="${iso}"></div><div id="l3Out" style="font-size:13.5px;color:#cfd8e6"></div><div class="checks"><label><input type="checkbox" id="s3Poly" checked> Aristas del poliedro</label></div>`;
    panel.querySelectorAll('#l3G button').forEach(b=>b.onclick=()=>{ g=b.dataset.k; mode='d'; idx=0; Acc=SYM.I3(); sample(); buildUI(); draw(); });
    panel.querySelectorAll('#l3F button').forEach(b=>b.onclick=()=>{ mode=b.dataset.m; idx=+b.dataset.i; Acc=SYM.I3(); sample(); buildUI(); draw(); });
    panel.querySelectorAll('#l3O button').forEach(b=>b.onclick=()=>{ if(b.dataset.o==='r'){ Acc=SYM.I3(); draw(); document.querySelector('#l3Out').innerHTML=''; return; } const op=G[g].ops[+b.dataset.o]; const r=SYM.analyze(g,mode,idx,Acc,op); const t0=performance.now();
      const step=now=>{ const t=Math.min(1,(now-t0)/1100); anim={op,t}; draw(); if(t<1) requestAnimationFrame(step); else { Acc=SYM.mul(SYM.finalMatrix(op),Acc); anim=null; draw(); document.querySelector('#l3Out').innerHTML=`<p><b>${op.label}</b>: la función ${r.kind}. χ = ${String(r.chi).replace('.',',')}; carácter de la representación ${r.irrep} en esta clase: <b>${r.table}</b>.</p>`; } }; requestAnimationFrame(step); });
    document.querySelector('#l3I').oninput=e=>{ iso=+e.target.value; document.querySelector('#l3IV').textContent=iso; sample(); draw(); }; document.querySelector('#s3Poly').onchange=draw; }
  cv.addEventListener('pointerdown',e=>{ drag=[e.clientX,e.clientY]; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener('pointermove',e=>{ if(!drag) return; yaw+=(e.clientX-drag[0])*0.01; pitch=Math.max(-1.4,Math.min(1.4,pitch+(e.clientY-drag[1])*0.01)); drag=[e.clientX,e.clientY]; draw(); });
  cv.addEventListener('pointerup',()=>drag=null); window.addEventListener('resize',draw);
  const el=document.getElementById('s3EList'); if(el) el.innerHTML='<p class="s3k">En modo sin conexión los elementos de simetría se aplican con los botones de operaciones del panel.</p>';
  sample(); buildUI(); draw(); }
