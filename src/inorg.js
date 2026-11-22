/* ===== Inorgánica I–III y temas emergentes: datos y cálculos (sin DOM) ===== */
const INORG=(function(){
"use strict";
/* ---------- isomería: enumeración por grupos de rotación ---------- */
function closure(gens,n){ const id=[...Array(n).keys()], key=p=>p.join(','), seen=new Map([[key(id),id]]), q=[id];
  while(q.length){ const p=q.pop(); gens.forEach(g=>{ const c=p.map(i=>g[i]); const k=key(c); if(!seen.has(k)){ seen.set(k,c); q.push(c); } }); } return [...seen.values()]; }
// octaedro: 0 +x, 1 −x, 2 +y, 3 −y, 4 +z, 5 −z
const OCT={n:6,pos:[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],rot:closure([[2,3,1,0,4,5],[0,1,4,5,3,2]],6),inv:[1,0,3,2,5,4],opp:[1,0,3,2,5,4],name:'octaédrica'};
// plano cuadrado: 0 +x, 1 +y, 2 −x, 3 −y (rotaciones propias en 3D, incluidas las C2 en el plano)
const SQP={n:4,pos:[[1,0,0],[0,1,0],[-1,0,0],[0,-1,0]],rot:closure([[1,2,3,0],[0,3,2,1]],4),inv:[0,3,2,1],opp:[2,3,0,1],name:'plano cuadrada',planar:true};
// tetraedro
const s3=1/Math.sqrt(3);
const TET={n:4,pos:[[s3,s3,s3],[-s3,-s3,s3],[-s3,s3,-s3],[s3,-s3,-s3]],rot:closure([[1,2,0,3],[0,2,3,1]],4),inv:[1,0,2,3],opp:[-1,-1,-1,-1],name:'tetraédrica'};
function parseIso(f){ let s=f.replace(/\s/g,'').replace(/^M/,''); const mono=[], bi=[]; let m;
  const re=/\(AA\)(\d*)|([a-z])(\d*)/g; while((m=re.exec(s))){ if(m[0]==='') break; if(m[0].startsWith('(')) { const k=m[1]?+m[1]:1; for(let i=0;i<k;i++) bi.push('A'); } else { const k=m[3]?+m[3]:1; for(let i=0;i<k;i++) mono.push(m[2]); } }
  return {mono,bi,cn:mono.length+2*bi.length}; }
function canon(arr,G){ let best=null; G.rot.forEach(g=>{ const b=new Array(G.n); arr.forEach((x,i)=>b[g[i]]=x); const map={}; let c=0; const s=b.map(x=>{ if(x[0]!=='#') return x; if(!(x in map)) map[x]='#'+(c++); return map[x]; }).join('|'); if(best===null||s<best) best=s; }); return best; }
function enumerate(f,G){ const {mono,bi,cn}=parseIso(f); if(cn!==G.n) return null;
  const found=new Map(); const arr=new Array(G.n).fill(null);
  const placeMono=()=>{ const free=[...arr.keys()].filter(i=>arr[i]===null); const ms=mono.slice().sort();
    const rec=(k,used)=>{ if(k===free.length){ const key=canon(arr,G); if(!found.has(key)) found.set(key,arr.slice()); return; } const seenL=new Set(); ms.forEach((l,j)=>{ if(used[j]||seenL.has(l)) return; seenL.add(l); used[j]=true; arr[free[k]]=l; rec(k+1,used); used[j]=false; arr[free[k]]=null; }); };
    rec(0,ms.map(()=>false)); };
  const placeBi=k=>{ if(k===bi.length){ placeMono(); return; } for(let i=0;i<G.n;i++) for(let j=i+1;j<G.n;j++){ if(arr[i]!==null||arr[j]!==null) continue; if(G.opp[i]===j) continue; if(G===TET) {} arr[i]=arr[j]='#'+k; placeBi(k+1); arr[i]=arr[j]=null; } };
  placeBi(0);
  const list=[...found.entries()].map(([key,a])=>{ const m=new Array(G.n); a.forEach((x,i)=>m[G.inv[i]]=x); const mk=canon(m,G); return {key,a,chiral:!G.planar&&mk!==key,mirror:mk}; });
  // agrupar enantiómeros
  const out=[], done=new Set(); list.forEach(o=>{ if(done.has(o.key)) return; done.add(o.key); if(o.chiral){ const p=list.find(x=>x.key===o.mirror); if(p) done.add(p.key); out.push({a:o.a,pair:p?p.a:null,chiral:true}); } else out.push({a:o.a,chiral:false}); });
  return {geo:out,total:list.length,cn};
}
function descIso(a,G){ // descripción cis/trans, fac/mer
  if(G!==OCT&&G!==SQP) return ''; const labs=[...new Set(a.filter(x=>x[0]!=='#'))]; const d=[];
  labs.forEach(l=>{ const idx=a.map((x,i)=>x===l?i:-1).filter(i=>i>=0); if(idx.length===2) d.push(`${l}: ${G.opp[idx[0]]===idx[1]?'trans':'cis'}`); if(idx.length===3&&G===OCT){ const tr=idx.some(i=>idx.includes(G.opp[i])); d.push(`${l}₃: ${tr?'mer':'fac'}`); } });
  return d.join(' · '); }

/* ---------- Latimer y Frost (medio ácido, V) ---------- */
const LAT={Mn:{sp:[['MnO₄⁻',7],['MnO₄²⁻',6],['MnO₂',4],['Mn³⁺',3],['Mn²⁺',2],['Mn',0]],E:[0.56,2.26,0.95,1.51,-1.18],el:'Mn'},
  Cl:{sp:[['ClO₄⁻',7],['ClO₃⁻',5],['HClO₂',3],['HClO',1],['Cl₂',0],['Cl⁻',-1]],E:[1.20,1.18,1.67,1.63,1.36],el:'Cl₂'},
  Fe:{sp:[['FeO₄²⁻',6],['Fe³⁺',3],['Fe²⁺',2],['Fe',0]],E:[2.20,0.77,-0.44],el:'Fe'},
  Cu:{sp:[['Cu²⁺',2],['Cu⁺',1],['Cu',0]],E:[0.15,0.52],el:'Cu'},
  O:{sp:[['O₂',0],['H₂O₂',-1],['H₂O',-2]],E:[0.70,1.76],el:'O₂'}};
function frost(key){ const L=LAT[key], sp=L.sp, i0=sp.findIndex(s=>s[1]===0); const Y=new Array(sp.length).fill(0);
  for(let i=i0-1;i>=0;i--){ const n=sp[i][1]-sp[i+1][1]; Y[i]=Y[i+1]+n*L.E[i]; }
  for(let i=i0+1;i<sp.length;i++){ const n=sp[i-1][1]-sp[i][1]; Y[i]=Y[i-1]-n*L.E[i-1]; }
  const dis=sp.map((s,i)=>i>0&&i<sp.length-1&&L.E[i]>L.E[i-1]);
  return {pts:sp.map((s,i)=>({sp:s[0],N:s[1],nE:Y[i],dis:dis[i]})),L}; }
function skipE(key,i,j){ // E° entre especies no adyacentes (i más oxidada)
  const L=LAT[key]; let num=0; for(let k=i;k<j;k++) num+=(L.sp[k][1]-L.sp[k+1][1])*L.E[k]; return num/(L.sp[i][1]-L.sp[j][1]); }

/* ---------- Pourbaix ---------- */
const S59=0.05916;
function pourbaixFe(){ // [Fe] = 1e-6 M; Ksp Fe(OH)3 2.8e-39, Fe(OH)2 4.9e-17
  const pH3=14+Math.log10(Math.cbrt(2.8e-39/1e-6)), pH2=14+Math.log10(Math.sqrt(4.9e-17/1e-6)), Efe=-0.44+S59/2*Math.log10(1e-6), E32=0.77;
  const E_at=pH=>E32-3*S59*(pH-pH3); const Eb=E_at(pH2);
  return {pH3,pH2,lines:[[[0,E32],[pH3,E32]],[[pH3,E32],[pH2,Eb]],[[pH2,Eb],[14,Eb-S59*(14-pH2)]],[[0,Efe],[pH2,Efe]],[[pH2,Efe],[14,Efe-S59*(14-pH2)]],[[pH3,E32],[pH3,1.6]]],
    labels:[['Fe³⁺',1.4,1.25],['Fe²⁺',4,0.2],['Fe(OH)₃(s)',9,0.5],['Fe(OH)₂(s)',11.5,-0.55],['Fe(s)',5,-0.9]],note:'[Fe] = 10⁻⁶ M; K<sub>ps</sub> Fe(OH)₃ = 2,8 × 10⁻³⁹ y Fe(OH)₂ = 4,9 × 10⁻¹⁷; E°(Fe³⁺/Fe²⁺) = 0,77 V y E°(Fe²⁺/Fe) = −0,44 V. Diagrama simplificado (sin Fe₃O₄ ni hidroxocomplejos).'}; }
const AS5=[2.20,6.97,11.53], AS3=[9.22,12.13];
function alphas(pKs,pH){ const h=Math.pow(10,-pH), K=pKs.map(p=>Math.pow(10,-p)); const terms=[1]; let prod=1; K.forEach(k=>{ prod*=k; terms.push(prod/Math.pow(h,terms.length)); }); const s=terms.reduce((a,b)=>a+b,0); return terms.map(t=>t/s); }
function asBoundary(pH){ const a5=alphas(AS5,pH)[0], a3=alphas(AS3,pH)[0]; return 0.560+S59/2*Math.log10(a5/a3)-S59*pH; }
/* ---------- Ellingham (por mol de O₂; ΔH en kJ, ΔS en J/K) ---------- */
const ELL=[['2Ag₂O',-62.2,-133,'4Ag + O₂ → 2Ag₂O'],['2Cu₂O',-337.2,-152,'4Cu + O₂ → 2Cu₂O'],['2PbO',-438,-202,'2Pb + O₂ → 2PbO'],['2NiO',-479.4,-189,'2Ni + O₂ → 2NiO'],['2H₂O(g)',-483.6,-89,'2H₂ + O₂ → 2H₂O(g)'],['2FeO',-544,-138,'2Fe + O₂ → 2FeO'],['SnO₂',-577.6,-207,'Sn + O₂ → SnO₂'],['2ZnO',-701,-201,'2Zn + O₂ → 2ZnO'],['⅔Cr₂O₃',-759.8,-183,'4/3Cr + O₂ → 2/3Cr₂O₃'],['TiO₂',-944,-185,'Ti + O₂ → TiO₂'],['⅔Al₂O₃',-1117,-209,'4/3Al + O₂ → 2/3Al₂O₃'],['2MgO',-1203,-217,'2Mg + O₂ → 2MgO'],['CO₂',-393.5,3,'C + O₂ → CO₂'],['2CO',-221,179,'2C + O₂ → 2CO'],['2CO₂ (desde CO)',-566,-173,'2CO + O₂ → 2CO₂']];
/* ---------- Born–Haber (kJ/mol) ---------- */
const BH={NaCl:{f:'NaCl',dfH:-411.2,steps:[['Sublimación de Na',107.3],['1.ª ionización de Na',495.8],['½ disociación de Cl₂',121.7],['Afinidad electrónica de Cl',-348.6]],lit:787},
  KCl:{f:'KCl',dfH:-436.5,steps:[['Sublimación de K',89.0],['1.ª ionización de K',418.8],['½ disociación de Cl₂',121.7],['Afinidad electrónica de Cl',-348.6]],lit:717},
  LiF:{f:'LiF',dfH:-616.0,steps:[['Sublimación de Li',159.4],['1.ª ionización de Li',520.2],['½ disociación de F₂',79.4],['Afinidad electrónica de F',-328.2]],lit:1037},
  MgO:{f:'MgO',dfH:-601.6,steps:[['Sublimación de Mg',147.1],['1.ª + 2.ª ionización de Mg',2188.4],['½ disociación de O₂',249.2],['1.ª afinidad electrónica de O',-141.0],['2.ª afinidad electrónica de O (endotérmica)',798]],lit:3850},
  CaF2:{f:'CaF₂',dfH:-1228.0,steps:[['Sublimación de Ca',178.2],['1.ª + 2.ª ionización de Ca',1735.2],['Disociación de F₂',158.8],['2 × afinidad electrónica de F',-656.4]],lit:2630}};
const MAD=[['NaCl (sal gema)',1.74756,6],['CsCl',1.76267,8],['ZnS blenda',1.63806,4],['ZnS wurtzita',1.64132,4],['CaF₂ fluorita',2.51939,8],['TiO₂ rutilo',2.408,6]];
const BORN={He:5,Ne:7,Ar:9,Kr:10,Xe:12};
/* ---------- ácidos ---------- */
const OXO=[['HClO','Cl(OH)',0,1,7.53],['HClO₂','ClO(OH)',1,1,1.96],['HClO₃','ClO₂(OH)',2,1,-1],['HClO₄','ClO₃(OH)',3,1,-10],['HNO₂','NO(OH)',1,1,3.3],['HNO₃','NO₂(OH)',2,1,-1.4],['H₂SO₃','SO(OH)₂',1,2,1.8],['H₂SO₄','SO₂(OH)₂',2,2,-3],['H₃PO₄','PO(OH)₃',1,3,2.15],['H₃PO₃','HPO(OH)₂ (un H unido a P)',1,2,1.3],['H₃PO₂','H₂PO(OH) (dos H unidos a P)',1,1,1.2],['H₅IO₆','IO(OH)₅',1,5,3.3],['H₃BO₃','B(OH)₃ (ácido de Lewis)',0,3,9.24],['H₄SiO₄','Si(OH)₄',0,4,9.8],['H₂CO₃','CO(OH)₂',1,2,6.35]];
const HYDR=[['HF',3.2],['HCl',-7],['HBr',-9],['HI',-10],['H₂O',14.0],['H₂S',7.0],['H₂Se',3.9],['H₂Te',2.6]];
const AQUO=[['Fe³⁺',2.2],['Cr³⁺',4.0],['Al³⁺',5.0],['Cu²⁺',7.5],['Zn²⁺',9.0],['Fe²⁺',9.5],['Mg²⁺',11.4],['Ca²⁺',12.7],['Na⁺',14.2]];
const HSAB={duro:['H⁺','Li⁺','Na⁺','K⁺','Be²⁺','Mg²⁺','Ca²⁺','Al³⁺','Cr³⁺','Fe³⁺','Ti⁴⁺','BF₃','SO₃','CO₂','F⁻','OH⁻','H₂O','O²⁻','NH₃','Cl⁻','NO₃⁻','SO₄²⁻','CO₃²⁻','PO₄³⁻'],intermedio:['Fe²⁺','Co²⁺','Ni²⁺','Cu²⁺','Zn²⁺','Pb²⁺','Sn²⁺','SO₂','Br⁻','NO₂⁻','N₃⁻','piridina','SO₃²⁻'],blando:['Cu⁺','Ag⁺','Au⁺','Hg²⁺','Cd²⁺','Pt²⁺','Pd²⁺','Tl⁺','BH₃','I₂','I⁻','S²⁻','SCN⁻ (por S)','CN⁻','CO','PR₃','H⁻','RS⁻','C₂H₄']};
/* ---------- especiación ---------- */
const SPEC={'Cu²⁺–NH₃':{M:'Cu',L:'NH₃',logK:[4.04,3.43,2.80,1.48]},'Ag⁺–NH₃':{M:'Ag',L:'NH₃',logK:[3.31,3.91]},'Ni²⁺–NH₃':{M:'Ni',L:'NH₃',logK:[2.80,2.24,1.73,1.19,0.75,0.03]},'Ni²⁺–en':{M:'Ni',L:'en',logK:[7.52,6.32,4.49]}};
function speciation(logK,pL){ const L=Math.pow(10,-pL); let b=1; const t=[1]; logK.forEach(k=>{ b*=Math.pow(10,k); t.push(b*Math.pow(L,t.length)); }); const s=t.reduce((a,c)=>a+c,0); return t.map(x=>x/s); }
/* ---------- masas atómicas para economía atómica ---------- */
const MM={H:1.008,He:4.0026,Li:6.94,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Na:22.990,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,K:39.098,Ca:40.078,Ti:47.867,V:50.942,Cr:51.996,Mn:54.938,Fe:55.845,Co:58.933,Ni:58.693,Cu:63.546,Zn:65.38,Br:79.904,Ag:107.868,Sn:118.71,I:126.904,Au:196.967,Pb:207.2,W:183.84,Sb:121.76};
function molarMass(el){ let m=0; for(const [e,n] of Object.entries(el)){ if(!(e in MM)) return NaN; m+=MM[e]*n; } return m; }
/* ---------- semiconductores, Ln, perovskitas, baterías ---------- */
const SEMI=[['Si',1.12,'indirecto; base de la microelectrónica y de las celdas solares'],['Ge',0.66,'primeros transistores (1947)'],['GaAs',1.42,'directo; láseres y celdas de alta eficiencia'],['GaN',3.4,'LED azul (Nobel de Física 2014)'],['CdS',2.42,'amarillo de cadmio'],['HgS (cinabrio)',2.0,'bermellón: rojo'],['Cu₂O',2.1,'rojo; fotocatálisis'],['CdSe',1.74,'puntos cuánticos'],['TiO₂ (anatasa)',3.2,'blanco; fotocatálisis con UV'],['ZnO',3.37,'blanco; protector solar'],['PbS',0.41,'detectores de IR'],['Diamante',5.47,'aislante, incoloro']];
const LN=[['La',1.032],['Ce',1.01],['Pr',0.99],['Nd',0.983],['Pm',0.97],['Sm',0.958],['Eu',0.947],['Gd',0.938],['Tb',0.923],['Dy',0.912],['Ho',0.901],['Er',0.89],['Tm',0.88],['Yb',0.868],['Lu',0.861]];
const PEROV=[['SrTiO₃',1.44,0.605,1.40,'cúbica a temperatura ambiente'],['BaTiO₃',1.61,0.605,1.40,'t > 1: Ti «suelto» en su octaedro → tetragonal y ferroeléctrica'],['CaTiO₃',1.34,0.605,1.40,'ortorrómbica (octaedros inclinados); el mineral perovskita'],['LaAlO₃',1.36,0.535,1.40,'romboédrica, casi cúbica; sustrato de películas delgadas'],['KNbO₃',1.64,0.64,1.40,'ferroeléctrica'],['CsPbI₃',1.88,1.19,2.20,'t ≈ 0,85: la fase negra cúbica solo es estable en caliente; a temperatura ambiente pasa a una fase amarilla no perovskita'],['CH₃NH₃PbI₃ (MAPbI₃)',2.17,1.19,2.20,'perovskita híbrida de las celdas solares; radio efectivo del catión orgánico'],['HC(NH₂)₂PbI₃ (FAPbI₃)',2.53,1.19,2.20,'catión formamidinio: t más cercano a 1'],['CsPbBr₃',1.88,1.19,1.96,'nanocristales luminiscentes verdes']];
const BATT=[['LiCoO₂ (cátodo)',97.87,1,3.9,'solo ≈ 0,5 Li se extrae de forma reversible: ≈ 140 mAh/g prácticos'],['LiFePO₄ (cátodo)',157.76,1,3.4,'olivino; muy estable y seguro'],['LiMn₂O₄ (cátodo)',180.81,1,4.0,'espinela; barato, con disolución de Mn'],['LiNi⅓Mn⅓Co⅓O₂ (NMC111)',96.46,1,3.7,'laminar; vehículos eléctricos'],['Grafito → LiC₆ (ánodo)',72.07,1,0.1,'intercalación entre capas de grafeno'],['Li₄Ti₅O₁₂ (ánodo)',459.09,3,1.55,'«cero deformación», carga rápida'],['Silicio → Li₁₅Si₄ (ánodo)',28.085,3.75,0.4,'gran capacidad pero se expande ≈ 300 %'],['Litio metálico (ánodo)',6.94,1,0,'máxima capacidad; dendritas']];
return {OCT,SQP,TET,parseIso,enumerate,descIso,canon,LAT,frost,skipE,pourbaixFe,AS5,AS3,alphas,asBoundary,ELL,BH,MAD,BORN,OXO,HYDR,AQUO,HSAB,SPEC,speciation,MM,molarMass,SEMI,LN,PEROV,BATT,S59};
})();
if(typeof module!=='undefined') module.exports=INORG;
