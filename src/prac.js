/* ===== Motores para entrenadores, redox, termodinámica y periódica (sin DOM) ===== */
const PRAC = (function(){
"use strict";
/* ---------- indexación de difractogramas cúbicos ---------- */
const NLIST={P:[1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21,22,24,25,26,27,29,30,32],I:[2,4,6,8,10,12,14,16,18,20,22,24,26,30,32],F:[3,4,8,11,12,16,19,20,24,27,32,35,36,40]};
function hklOf(N,type){ const out=[]; for(let h=0;h*h<=N;h++)for(let k=0;k<=h;k++)for(let l=0;l<=k;l++){ if(h*h+k*k+l*l!==N) continue;
  if(type==='I'&&(h+k+l)%2) continue; if(type==='F'&&!((h%2===k%2)&&(k%2===l%2))) continue; out.push(`${h}${k}${l}`); } return out; }
function indexCubic(tt,lam){ const s2=tt.map(t=>Math.sin(t*Math.PI/360)**2).sort((a,b)=>a-b); const res=[];
  ['P','I','F'].forEach(type=>{ const L=NLIST[type]; for(let k0=0;k0<(type==='F'?2:1);k0++){ const s=s2[0]/L[k0]; let err=0, ok=true; const Ns=s2.map(x=>{ const n=x/s; const best=L.reduce((b,v)=>Math.abs(v-n)<Math.abs(b-n)?v:b,L[0]); err+=Math.abs(best-n)/n; return {n,N:best}; });
    // los N deben ser crecientes y sin saltarse demasiados valores permitidos
    const used=new Set(Ns.map(o=>o.N)); if(used.size<Ns.length) ok=false; const skipped=L.filter(v=>v<=Math.max(...used)&&!used.has(v)).length;
    const a=lam/(2*Math.sqrt(s)); res.push({type,k0,err:err/s2.length,skipped,ok,a,Ns:Ns.map((o,i)=>({tt:2*Math.asin(Math.sqrt(s2[i]))*180/Math.PI,s2:s2[i],n:o.n,N:o.N,hkl:hklOf(o.N,type)}))}); } });
  res.sort((x,y)=>(x.ok===y.ok?0:x.ok?-1:1)||(x.err+0.02*x.skipped)-(y.err+0.02*y.skipped)); return res; }
// celdas ortogonales: 1/d² = h²A + k²B + l²C (mínimos cuadrados)
function fitOrtho(rows,lam){ // rows: [h,k,l,2θ]
  const M=[[0,0,0],[0,0,0],[0,0,0]], v=[0,0,0];
  rows.forEach(([h,k,l,t])=>{ const d=lam/(2*Math.sin(t*Math.PI/360)), y=1/(d*d), x=[h*h,k*k,l*l]; for(let i=0;i<3;i++){ v[i]+=x[i]*y; for(let j=0;j<3;j++) M[i][j]+=x[i]*x[j]; } });
  const used=[0,1,2].filter(i=>M[i][i]>0); const sub=used.map(i=>used.map(j=>M[i][j])), b=used.map(i=>v[i]); const sol=solve(sub,b); const out=[null,null,null]; used.forEach((i,k)=>out[i]=sol?1/Math.sqrt(sol[k]):null); return out; }
function solve(A,b){ const n=A.length, M=A.map((r,i)=>r.concat([b[i]])); for(let c=0;c<n;c++){ let p=c; for(let r=c+1;r<n;r++) if(Math.abs(M[r][c])>Math.abs(M[p][c])) p=r; if(Math.abs(M[p][c])<1e-14) return null; [M[c],M[p]]=[M[p],M[c]];
  for(let r=0;r<n;r++){ if(r===c) continue; const f=M[r][c]/M[c][c]; for(let k=c;k<=n;k++) M[r][k]-=f*M[c][k]; } } return M.map((r,i)=>r[n]/r[i]); }

/* ---------- igualación de ecuaciones (método algebraico) ---------- */
function parseSpecies(str){ // "Na[Au(CN)2]", "IO3-", "Fe3+", "SO4 2-", "Cr2O7^2-", "SO₄²⁻"
  let s=str.trim().replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]+[⁺⁻]|[⁺⁻]/g,m=>'^'+m.replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]/g,c=>'⁰¹²³⁴⁵⁶⁷⁸⁹'.indexOf(c)).replace('⁺','+').replace('⁻','-'))
     .replace(/[₀-₉]/g,c=>String(c.charCodeAt(0)-8320)); let q=0;
  let m=s.match(/[\^\s]\(?(\d*)([+-])\)?$/);
  if(m){ q=(m[1]?+m[1]:1)*(m[2]==='+'?1:-1); s=s.slice(0,m.index).trim(); }
  else if((m=s.match(/(\d*)([+-])$/))){ const head=s.slice(0,m.index); const single=/^[A-Z][a-z]?$/.test(head);
    if(single&&m[1]){ q=+m[1]*(m[2]==='+'?1:-1); s=head; } else { q=m[2]==='+'?1:-1; s=s.slice(0,s.length-1); } }
  s=s.replace(/\[/g,'(').replace(/\]/g,')').replace(/·/g,'.');
  const parts=s.split('.'); const tot={};
  parts.forEach(p=>{ const mm=p.match(/^(\d+)(.*)$/); const mult=mm?+mm[1]:1; const body=mm?mm[2]:p; const c=parseGroup(body); Object.entries(c).forEach(([e,n])=>tot[e]=(tot[e]||0)+n*mult); });
  return {el:tot,q}; }
function parseGroup(s){ let i=0; const rec=()=>{ const out={}; while(i<s.length){ const ch=s[i];
    if(ch==='('){ i++; const inner=rec(); let n=''; while(/\d/.test(s[i]||'')) n+=s[i++]; const k=n?+n:1; Object.entries(inner).forEach(([e,v])=>out[e]=(out[e]||0)+v*k); }
    else if(ch===')'){ i++; return out; }
    else if(/[A-Z]/.test(ch)){ let e=ch; i++; while(/[a-z]/.test(s[i]||'')) e+=s[i++]; let n=''; while(/\d/.test(s[i]||'')) n+=s[i++]; out[e]=(out[e]||0)+(n?+n:1); }
    else i++; } return out; }; return rec(); }
function gcd(a,b){ a=Math.abs(a); b=Math.abs(b); while(b){ [a,b]=[b,a%b]; } return a; }
function balance(eq){ // "A + B -> C + D"
  const [L,R]=eq.split(/->|→|=/).map(x=>x.split(/\s\+\s/).map(t=>t.trim()).filter(Boolean)); if(!L||!R) return null;
  const sp=L.concat(R).map(parseSpecies); const els=[...new Set(sp.flatMap(s=>Object.keys(s.el)))]; const rows=els.map(e=>sp.map((s,j)=>(s.el[e]||0)*(j<L.length?1:-1)));
  if(sp.some(s=>s.q)) rows.push(sp.map((s,j)=>s.q*(j<L.length?1:-1)));
  // espacio nulo con fracciones
  const n=sp.length, F=(a,b=1)=>{ const g=gcd(a,b)||1; return b<0?[-a/g,-b/g]:[a/g,b/g]; }, mul=(x,y)=>F(x[0]*y[0],x[1]*y[1]), sub=(x,y)=>F(x[0]*y[1]-y[0]*x[1],x[1]*y[1]), div=(x,y)=>F(x[0]*y[1],x[1]*y[0]);
  const M=rows.map(r=>r.map(v=>[v,1])); let piv=[], r=0;
  for(let c=0;c<n&&r<M.length;c++){ let p=-1; for(let i=r;i<M.length;i++) if(M[i][c][0]!==0){ p=i; break; } if(p<0) continue; [M[r],M[p]]=[M[p],M[r]]; const pv=M[r][c]; M[r]=M[r].map(x=>div(x,pv));
    for(let i=0;i<M.length;i++){ if(i===r||M[i][c][0]===0) continue; const f=M[i][c]; M[i]=M[i].map((x,k)=>sub(x,mul(f,M[r][k]))); } piv.push(c); r++; }
  const free=[...Array(n).keys()].filter(c=>!piv.includes(c)); if(free.length!==1) return {err:free.length?'La ecuación admite varias combinaciones independientes: separa las reacciones o fija un coeficiente.':'No se puede igualar con estas especies.'};
  const f=free[0]; const x=new Array(n).fill(null); x[f]=[1,1]; piv.forEach((c,i)=>{ x[c]=F(-M[i][f][0],M[i][f][1]); });
  let L2=x.reduce((s,v)=>s*v[1]/gcd(s,v[1]),1); let coef=x.map(v=>v[0]*L2/v[1]); if(coef.every(c=>c<=0)) coef=coef.map(c=>-c); if(coef.some(c=>c<=0)) return {err:'Alguna especie queda con coeficiente cero o negativo: revisa que estén en el lado correcto.'};
  const g=coef.reduce((a,b)=>gcd(a,b)); coef=coef.map(c=>c/g); return {L,R,coef,sp}; }
// estados de oxidación por reglas
const FIX={F:-1,Li:1,Na:1,K:1,Rb:1,Cs:1,Be:2,Mg:2,Ca:2,Sr:2,Ba:2,Al:3,Zn:2,Ag:1};
function oxStates(str){ const {el,q}=parseSpecies(str); const els=Object.keys(el); const ox={};
  if(els.length===1){ ox[els[0]]=q/el[els[0]]; return ox; }
  let rest=q, unk=[]; const has=e=>el[e]>0;
  els.forEach(e=>{ if(FIX[e]!=null){ ox[e]=FIX[e]; rest-=FIX[e]*el[e]; } });
  if(has('H')&&ox.H==null){ const metalH=els.every(e=>e==='H'||['Li','Na','K','Ca','Mg','Al','B'].includes(e))&&els.length>1&&!has('O'); ox.H=metalH&&!has('B')?-1:1; rest-=ox.H*el.H; }
  if(has('O')&&ox.O==null){ ox.O=has('F')?null:-2; if(ox.O!=null) rest-=ox.O*el.O; }
  if(/CN/.test(str)&&has('C')&&has('N')&&ox.C==null){ ox.C=2; ox.N=-3; rest-=2*el.C-3*el.N; }
  const halo=['Cl','Br','I'].filter(e=>has(e)&&ox[e]==null);
  unk=els.filter(e=>ox[e]==null);
  if(unk.length>1&&halo.length){ halo.forEach(e=>{ if(unk.length>1&&!has('O')){ ox[e]=-1; rest-=-1*el[e]; unk=unk.filter(x=>x!==e); } }); }
  if(unk.length===1){ const e=unk[0]; ox[e]=rest/el[e]; } else unk.forEach(e=>ox[e]=NaN);
  return ox; }
function redoxInfo(b){ if(!b||b.err) return null; const n=b.L.length; const before={}, after={};
  b.L.forEach(s=>{ const o=oxStates(s); Object.entries(o).forEach(([e,v])=>{ (before[e]=before[e]||[]).push([s,v]); }); });
  b.R.forEach(s=>{ const o=oxStates(s); Object.entries(o).forEach(([e,v])=>{ (after[e]=after[e]||[]).push([s,v]); }); });
  const changes=[]; Object.keys(before).forEach(e=>{ before[e].forEach(([s,v])=>{ (after[e]||[]).forEach(([s2,v2])=>{ if(isFinite(v)&&isFinite(v2)&&Math.abs(v-v2)>1e-9) changes.push({el:e,from:s,to:s2,a:v,b:v2,type:v2>v?'oxidación':'reducción'}); }); }); });
  const uniq=[]; changes.forEach(c=>{ if(!uniq.some(u=>u.el===c.el&&u.from===c.from&&u.to===c.to)) uniq.push(c); });
  return {changes:uniq,redox:uniq.length>0,oxidant:[...new Set(uniq.filter(c=>c.type==='reducción').map(c=>c.from))],reductant:[...new Set(uniq.filter(c=>c.type==='oxidación').map(c=>c.from))]}; }

/* ---------- funciones de onda hidrogenoides (forma del curso) ---------- */
const WF=[
 {f:'e<sup>−Zr/a₀</sup>',o:'1s',n:1,l:0},{f:'(2 − Zr/a₀)·e<sup>−Zr/2a₀</sup>',o:'2s',n:2,l:0},{f:'r·e<sup>−Zr/2a₀</sup>·cos θ',o:'2p<sub>z</sub>',n:2,l:1},{f:'r·e<sup>−Zr/2a₀</sup>·sen θ·cos φ',o:'2p<sub>x</sub>',n:2,l:1},{f:'r·e<sup>−Zr/2a₀</sup>·sen θ·sen φ',o:'2p<sub>y</sub>',n:2,l:1},
 {f:'(27 − 18Zr/a₀ + 2Z²r²/a₀²)·e<sup>−Zr/3a₀</sup>',o:'3s',n:3,l:0},{f:'(6 − Zr/a₀)·r·e<sup>−Zr/3a₀</sup>·cos θ',o:'3p<sub>z</sub>',n:3,l:1},{f:'r²·e<sup>−Zr/3a₀</sup>·(3cos²θ − 1)',o:'3d<sub>z²</sub>',n:3,l:2},{f:'r²·e<sup>−Zr/3a₀</sup>·sen θ·cos θ·cos φ',o:'3d<sub>xz</sub>',n:3,l:2},{f:'r²·e<sup>−Zr/3a₀</sup>·sen θ·cos θ·sen φ',o:'3d<sub>yz</sub>',n:3,l:2},
 {f:'r²·e<sup>−Zr/3a₀</sup>·sen²θ·cos 2φ',o:'3d<sub>x²−y²</sub>',n:3,l:2},{f:'r²·e<sup>−Zr/3a₀</sup>·sen²θ·sen 2φ',o:'3d<sub>xy</sub>',n:3,l:2},{f:'r·e<sup>−Zr/2a₀</sup>·sen θ·e<sup>iφ</sup>',o:'2p<sub>+1</sub> (función compleja, m = +1)',n:2,l:1,m:1},{f:'r²·e<sup>−Zr/3a₀</sup>·sen²θ·e<sup>−2iφ</sup>',o:'3d<sub>−2</sub> (función compleja, m = −2)',n:3,l:2,m:-2}];

/* ---------- datos periódicos (Z = 1–36) ---------- */
const SYM=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'];
const IE1=[13.598,24.587,5.392,9.323,8.298,11.260,14.534,13.618,17.423,21.565,5.139,7.646,5.986,8.152,10.487,10.360,12.968,15.760,4.341,6.113,6.561,6.828,6.746,6.767,7.434,7.902,7.881,7.640,7.726,9.394,5.999,7.899,9.789,9.752,11.814,14.000];
const ENP=[2.20,null,0.98,1.57,2.04,2.55,3.04,3.44,3.98,null,0.93,1.31,1.61,1.90,2.19,2.58,3.16,null,0.82,1.00,1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65,1.81,2.01,2.18,2.55,2.96,3.00];
const MP=[-259,-272,181,1287,2076,3550,-210,-218,-220,-249,98,650,660,1414,44,115,-101,-189,64,842];
/* ---------- termodinámica ---------- */
const HYD=[['Li⁺',1,0.76,-519],['Na⁺',1,1.02,-406],['K⁺',1,1.38,-322],['Rb⁺',1,1.52,-293],['Cs⁺',1,1.67,-264],['Mg²⁺',2,0.72,-1921],['Ca²⁺',2,1.00,-1577],['Sr²⁺',2,1.18,-1443],['Ba²⁺',2,1.35,-1305],['Al³⁺',3,0.535,-4660],['In³⁺',3,0.80,-4109]];
const EPOT=[['F₂ + 2e⁻ → 2F⁻',2.87,'F₂','F⁻',2],['Au³⁺ + 3e⁻ → Au',1.50,'Au³⁺','Au',3],['Cl₂ + 2e⁻ → 2Cl⁻',1.36,'Cl₂','Cl⁻',2],['O₂ + 4H⁺ + 4e⁻ → 2H₂O',1.23,'O₂ (ácido)','H₂O',4],['MnO₂ + 4H⁺ + 2e⁻ → Mn²⁺ + 2H₂O',1.23,'MnO₂','Mn²⁺',2],['Br₂ + 2e⁻ → 2Br⁻',1.07,'Br₂','Br⁻',2],['Ag⁺ + e⁻ → Ag',0.80,'Ag⁺','Ag',1],['Fe³⁺ + e⁻ → Fe²⁺',0.77,'Fe³⁺','Fe²⁺',1],['I₂ + 2e⁻ → 2I⁻',0.54,'I₂','I⁻',2],['Cu²⁺ + 2e⁻ → Cu',0.34,'Cu²⁺','Cu',2],['2H⁺ + 2e⁻ → H₂',0.00,'H⁺','H₂',2],['Pb²⁺ + 2e⁻ → Pb',-0.13,'Pb²⁺','Pb',2],['Fe²⁺ + 2e⁻ → Fe',-0.44,'Fe²⁺','Fe',2],['Zn²⁺ + 2e⁻ → Zn',-0.76,'Zn²⁺','Zn',2],['Al³⁺ + 3e⁻ → Al',-1.66,'Al³⁺','Al',3],['Mg²⁺ + 2e⁻ → Mg',-2.37,'Mg²⁺','Mg',2],['Na⁺ + e⁻ → Na',-2.71,'Na⁺','Na',1],['Li⁺ + e⁻ → Li',-3.04,'Li⁺','Li',1]];
const SIGNS=[['Evaporación del agua a 25 °C','+','+','+'],['Congelación del agua a −10 °C','−','−','−'],['N₂ + 3H₂ → 2NH₃ (25 °C)','−','−','−'],['Cl₂ → 2Cl•','+','+','+'],['2H• → H₂','−','−','−'],['Hidratación de un catión, M⁺(g) → M⁺(ac)','−','−','−'],['CaCO₃ → CaO + CO₂ a 25 °C','+','+','+'],['Disolución de NH₄NO₃ en agua (se enfría)','−','+','+'],['Combustión del metano','−','−','−'],['C(grafito) → C(diamante) a 25 °C','+','−','+']]; // ΔG, ΔS, ΔH
/* ---------- llama y emisión ---------- */
const FLAME=[['Li','rojo carmín',[[670.8,1]],'#d0133d'],['Na','amarillo',[[589.0,1],[589.6,0.5]],'#ffc400'],['K','lila (violeta pálido)',[[766.5,1],[769.9,0.6],[404.4,0.1]],'#b48be0'],['Ca','rojo anaranjado',[[422.7,0.5],[622,0.9],[554,0.4]],'#ff6a2a'],['Sr','rojo intenso',[[460.7,0.4],[606,0.6],[650,0.9],[682,0.7]],'#e0102a'],['Ba','verde manzana',[[553.5,1],[513,0.5],[524,0.5]],'#9ad52a'],['Cu','verde azulado',[[510.6,0.5],[515.3,0.6],[521.8,0.8]],'#20c0a0']];
return {NLIST,hklOf,indexCubic,fitOrtho,parseSpecies,balance,oxStates,redoxInfo,WF,SYM,IE1,ENP,MP,HYD,EPOT,SIGNS,FLAME,gcd};
})();
if(typeof module!=='undefined') module.exports=PRAC;
