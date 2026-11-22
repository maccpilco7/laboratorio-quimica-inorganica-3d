/* ===== Átomo, núcleo y espectroscopía: cálculos puros (sin DOM) ===== */
const ATOM = (function(){
"use strict";
const h=6.62607015e-34, c=2.99792458e8, NA=6.02214076e23, eV=1.602176634e-19, RH=1.097e7; // R usado en los exámenes del curso
/* ---------- espectro del hidrógeno y de iones hidrogenoides ---------- */
const SERIES={1:'Lyman',2:'Balmer',3:'Paschen',4:'Brackett',5:'Pfund',6:'Humphreys'};
function region(nm){ if(nm<10) return 'rayos X'; if(nm<380) return 'ultravioleta'; if(nm<=750) return 'visible'; if(nm<1e6) return 'infrarrojo'; return 'microondas'; }
function rydberg(Z,ni,nf,R=RH){ // ni > nf: emisión
  const inv=nf===Infinity?0:1/(nf*nf), ini=ni===Infinity?0:1/(ni*ni);
  const wn=R*Z*Z*Math.abs(inv-ini); const E=h*c*wn; const lam=1/wn;
  return {wn, wnCm:wn/100, E, eV:E/eV, kJ:E*NA/1000, lam, nm:lam*1e9, nu:c*lam**-1, series:SERIES[Math.min(ni,nf)]||('n = '+Math.min(ni,nf)), region:region(lam*1e9)};
}
function ionization(Z,n,R=RH){ return rydberg(Z,n,Infinity,R); }
/* ---------- términos espectroscópicos por recuento de microestados ---------- */
const LL='SPDFGHIKLMNOQ';
const LSYM=['s','p','d','f'];
function comb(n,k){ if(k<0||k>n) return 0; let r=1; for(let i=1;i<=k;i++) r=r*(n-k+i)/i; return Math.round(r); }
function shellStates(l,n){ // microestados de una subcapa l con n electrones: lista de {ML, MS2, occ}
  const orbs=[]; for(let m=l;m>=-l;m--){ orbs.push([m,1]); orbs.push([m,-1]); }
  const out=[]; const N=orbs.length;
  (function rec(start,left,ML,MS2,sel){ if(!left){ out.push({ML,MS2,sel:sel.slice()}); return; }
    for(let i=start;i<=N-left;i++){ sel.push(i); rec(i+1,left-1,ML+orbs[i][0],MS2+orbs[i][1],sel); sel.pop(); } })(0,n,0,0,[]);
  return out;
}
function termsFromShells(shells){ // shells: [{l,n}]
  let tab=new Map(); tab.set('0,0',1);
  shells.forEach(({l,n})=>{ const st=shellStates(l,n); const nt=new Map(); tab.forEach((cnt,key)=>{ const [a,b]=key.split(',').map(Number); st.forEach(s=>{ const k2=(a+s.ML)+','+(b+s.MS2); nt.set(k2,(nt.get(k2)||0)+cnt); }); }); tab=nt; });
  const total=[...tab.values()].reduce((s,x)=>s+x,0);
  const terms=[]; const t=new Map(tab);
  for(;;){ let best=null; t.forEach((cnt,key)=>{ if(cnt<=0) return; const [ML,MS2]=key.split(',').map(Number); if(!best||ML>best[0]||(ML===best[0]&&MS2>best[1])) best=[ML,MS2]; });
    if(!best) break; const [L,S2]=best; if(L<0||S2<0) break;
    for(let ml=-L;ml<=L;ml++) for(let ms=-S2;ms<=S2;ms+=2){ const k=ml+','+ms; t.set(k,t.get(k)-1); }
    const Js=[]; for(let J2=Math.abs(2*L-S2);J2<=2*L+S2;J2+=2) Js.push(J2);
    terms.push({L,S2,mult:S2+1,sym:(S2+1)+LL[L],Js,deg:(2*L+1)*(S2+1)}); }
  return {terms,total};
}
const J2s=J2=>J2%2?`${J2}/2`:String(J2/2);
function hundGround(shells){ // reglas de Hund; la regla de J se aplica a la subcapa abierta con l > 0 (o a la única abierta)
  let S2=0,L=0; const boxes=[];
  shells.forEach(({l,n})=>{ const k=2*l+1; const up=Math.min(n,k), dn=n-up; const ms=[]; for(let i=0;i<k;i++) ms.push(l-i);
    let ML=0; for(let i=0;i<up;i++) ML+=ms[i]; for(let i=0;i<dn;i++) ML+=ms[i]; L+=ML; S2+=up-dn; boxes.push({l,n,ms,up,dn}); });
  // J: la subcapa que decide es la de mayor l; si está más que semillena → J máximo
  const key=shells.slice().sort((a,b)=>b.l-a.l)[0]; const half=2*key.l+1;
  const more=key.n>half; const J2=more?2*L+S2:Math.abs(2*L-S2);
  return {L,S2,J2,sym:`${S2+1}${LL[L]}${J2s(J2)}`,boxes,rule:more?'subcapa más que semillena → J = L + S':'subcapa semillena o menos → J = |L − S|'};
}
function parseConfig(str){ // "3d5", "1s1 2s1", "2p4 3s1", "[Ar]3d7"; solo subcapas abiertas
  const re=/(\d)([spdf])(\d{1,2})/g; let m; const out=[]; while((m=re.exec(str))){ const l=LSYM.indexOf(m[2]), n=+m[3]; if(n>2*(2*l+1)) return null; out.push({n:+m[1],sub:m[2],l,e:n}); }
  return out; }
function landeG(L,S2,J2){ const S=S2/2,J=J2/2; if(J===0) return 0; return 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1)); }
function zeeman(up,lo){ // up/lo: {L,S2,J2}
  const gu=landeG(up.L,up.S2,up.J2), gl=landeG(lo.L,lo.S2,lo.J2), lines=[];
  for(let mu=-up.J2;mu<=up.J2;mu+=2) for(let ml=-lo.J2;ml<=lo.J2;ml+=2){ const d=(mu-ml)/2; if(Math.abs(d)<=1){ if(up.J2===lo.J2&&mu===0&&ml===0&&up.J2%2===0) continue; lines.push({MJu:mu,MJl:ml,dM:d,shift:(gu*mu-gl*ml)/2}); } }
  return {gu,gl,lines};
}
// multiplicidad del estado fundamental de un complejo octaédrico dⁿ
function octSpin(n){ const hs=[0,1,2,3,4,5,4,3,2,1,0][n], ls=[0,1,2,3,2,1,0,1,2,1,0][n]; return {hs:hs+1,ls:ls+1,unpHS:hs,unpLS:ls}; }
/* ---------- operadores sobre ψnlm ---------- */
function operators(n,l,m,Z=1){ const ok=n>=1&&l>=0&&l<n&&Math.abs(m)<=l;
  return {ok,E:-13.606*Z*Z/(n*n),L2:l*(l+1),Lz:m,orb:ok?`${n}${'spdfghik'[l]||'?'}`:null,radialNodes:n-l-1,angNodes:l,why:ok?'':(l>=n?`l debe ser menor que n (l ≤ ${n-1})`:`|m| debe ser ≤ l (= ${l})`)}; }
/* ---------- núcleo ---------- */
const AR={H:1.008,He:4.0026,Li:6.94,Be:9.012,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Na:22.99,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,Ar:39.95,K:39.098,Ca:40.078,Co:58.933,Br:79.904,I:126.90,Cs:132.91,Pb:207.2,Bi:208.98,Po:209,Rn:222,Ra:226,Th:232.04,U:238.03,Tc:98,Sr:87.62};
const ZOF={H:1,He:2,Li:3,Be:4,B:5,C:6,N:7,O:8,F:9,Na:11,Mg:12,Al:13,Si:14,P:15,S:16,Cl:17,Ar:18,K:19,Ca:20,Co:27,Br:35,Sr:38,Tc:43,I:53,Cs:55,Pb:82,Bi:83,Po:84,Rn:86,Ra:88,Th:90,U:92};
function decayGuess(el,A){ const Z=ZOF[el], ar=AR[el]; if(!Z) return null;
  if(Z>83) return {mode:'α (y a veces β⁻)',why:'Z > 83: ningún isótopo es estable; los núcleos pesados pierden masa emitiendo partículas α.'};
  const ref=Math.round(ar); if(A<ref) return {mode:'β⁺ o captura electrónica',why:`A = ${A} es menor que la masa atómica media (${ar}): exceso de protones (N/Z bajo). Un protón se convierte en neutrón: p → n + e⁺ + ν.`};
  if(A>ref) return {mode:'β⁻',why:`A = ${A} es mayor que la masa atómica media (${ar}): exceso de neutrones (N/Z alto). Un neutrón se convierte en protón: n → p + e⁻ + ν̄.`};
  return {mode:'probablemente estable',why:'A coincide con el isótopo más abundante.'};
}
function seriesCount(Zi,Ai,Zf,Af){ const a=(Ai-Af)/4; const b=2*a-(Zi-Zf); return {alpha:a,beta:b,ok:Number.isInteger(a)&&b>=0}; }
const YEAR=365.25*24*3600;
function decay(t12,t){ const k=Math.LN2/t12; return {k,frac:Math.exp(-k*t)}; }
function halfLifeFromActivity(A0,A,t){ return t*Math.LN2/Math.log(A0/A); }
function bindingEnergy(Z,A,Matom){ const mH=1.00782503, mn=1.00866492, uc2=931.494; const dm=Z*mH+(A-Z)*mn-Matom; const MeV=dm*uc2;
  const Jmol=dm*1e-3*c*c; return {dm,MeV,perNucleon:MeV/A,kJg:Jmol/Matom/1000,Jmol}; }
/* ---------- RMN ---------- */
const NUC={ '1H':{I:0.5,g:42.577}, '2H':{I:1,g:6.536}, '13C':{I:0.5,g:10.708}, '19F':{I:0.5,g:40.078}, '31P':{I:0.5,g:17.235}, '15N':{I:0.5,g:-4.316}, '35Cl':{I:1.5,g:4.176}, '11B':{I:1.5,g:13.663}, '17O':{I:2.5,g:-5.774}, '29Si':{I:0.5,g:-8.465} };
function spinStates(I){ const out=[]; for(let m2=2*I;m2>=-2*I;m2-=2){ const m=m2/2; out.push({m,theta:Math.acos(m/Math.sqrt(I*(I+1)))*180/Math.PI}); } return out; }
/* ---------- masas: fórmula → cúmulo isotópico ---------- */
const ISO={H:[[1.007825,99.9885],[2.014102,0.0115]],D:[[2.014102,100]],B:[[10.012937,19.9],[11.009305,80.1]],C:[[12,98.93],[13.003355,1.07]],N:[[14.003074,99.636],[15.000109,0.364]],
  O:[[15.994915,99.757],[16.999132,0.038],[17.99916,0.205]],F:[[18.998403,100]],Si:[[27.976927,92.223],[28.976495,4.685],[29.97377,3.092]],P:[[30.973762,100]],
  S:[[31.972071,94.99],[32.971459,0.75],[33.967867,4.25],[35.967081,0.01]],Cl:[[34.968853,75.76],[36.965903,24.24]],Br:[[78.918337,50.69],[80.916291,49.31]],I:[[126.904473,100]]};
const ISO_SIMPLE={Cl:[[35,75],[37,25]],Br:[[79,50],[81,50]],B:[[10,20],[11,80]],H:[[1,100]],C:[[12,100]],O:[[16,100]],N:[[14,100]],S:[[32,100]],F:[[19,100]],D:[[2,100]],P:[[31,100]],I:[[127,100]],Si:[[28,100]]};
function parseFormula(f){ const re=/([A-Z][a-z]?)(\d*)/g; let m; const out={}; let ok=true; const s=f.replace(/[^A-Za-z0-9]/g,'');
  let consumed=0; while((m=re.exec(s))){ if(!m[0]) break; consumed+=m[0].length; if(!ISO[m[1]]) ok=false; out[m[1]]=(out[m[1]]||0)+(m[2]?+m[2]:1); }
  return ok&&consumed===s.length&&Object.keys(out).length?out:null; }
function cluster(counts,simple){ let dist=new Map([[0,{p:1,mass:0}]]);
  Object.entries(counts).forEach(([el,k])=>{ const iso=(simple&&ISO_SIMPLE[el])||ISO[el]; for(let i=0;i<k;i++){ const nd=new Map(); dist.forEach((v,key)=>iso.forEach(([ms,ab])=>{ const nk=key+Math.round(ms); const e=nd.get(nk)||{p:0,mass:0}; const pp=v.p*ab/100; e.mass=(e.mass*e.p+(v.mass+ms)*pp)/(e.p+pp||1); e.p+=pp; nd.set(nk,e); })); dist=new Map([...nd].filter(([,v])=>v.p>1e-9)); } });
  const arr=[...dist].map(([k,v])=>({mz:k,p:v.p,mass:v.mass})).sort((a,b)=>a.mz-b.mz); const mx=Math.max(...arr.map(x=>x.p)); return arr.map(x=>({mz:x.mz,rel:100*x.p/mx,abs:x.p,mass:x.mass})).filter(x=>x.rel>=0.1); }
function fragments(counts,simple){ const out=[{name:'M⁺',counts}]; Object.keys(counts).forEach(el=>{ if(counts[el]>=1){ const c2=Object.assign({},counts); c2[el]--; if(!c2[el]) delete c2[el]; if(Object.keys(c2).length) out.push({name:`[M − ${el}]⁺`,counts:c2}); } });
  Object.keys(counts).filter(e=>e!=='H').forEach(el=>out.push({name:`${el}⁺`,counts:{[el]:1}}));
  return out.map(f=>({name:f.name,formula:Object.entries(f.counts).map(([e,k])=>e+(k>1?k:'')).join(''),peaks:cluster(f.counts,simple)})); }
/* ---------- OM de diatómicas del 2.º periodo y HF ---------- */
const VE={H:1,Li:1,Be:2,B:3,C:4,N:5,O:6,F:7,Ne:8};
function diatomic(A,B,q){ // q = carga
  if(A==='H'||B==='H'){ // HF, HCl simplificado
    const X=A==='H'?B:A; const ne=1+VE[X]-q;
    const lev=[{n:'2s (no enl.)',t:'nb',k:'n',cap:2},{n:'σ',t:'σ',k:'b',cap:2},{n:'2pₓ, 2pᵧ (no enl.)',t:'nb',k:'n',cap:4},{n:'σ*',t:'σ',k:'a',cap:2}];
    return fill(lev,ne,`${A}${B}`); }
  const ne=VE[A]+VE[B]-q; const mix=(VE[A]+VE[B])<=10; // Li2..N2: σ2p por encima de π2p
  const lev=[{n:'σ2s',t:'σ',k:'b',cap:2},{n:'σ*2s',t:'σ',k:'a',cap:2}].concat(mix?[{n:'π2p',t:'π',k:'b',cap:4},{n:'σ2p',t:'σ',k:'b',cap:2}]:[{n:'σ2p',t:'σ',k:'b',cap:2},{n:'π2p',t:'π',k:'b',cap:4}]).concat([{n:'π*2p',t:'π',k:'a',cap:4},{n:'σ*2p',t:'σ',k:'a',cap:2}]);
  return Object.assign(fill(lev,ne,A+B),{mix}); }
function fill(lev,ne,name){ let left=ne; let unp=0; lev.forEach(L=>{ L.e=Math.min(L.cap,left); left-=L.e; const orbs=L.cap/2; L.unp=L.e<=orbs?L.e:L.cap-L.e; unp+=L.unp; });
  const b=lev.filter(L=>L.k==='b').reduce((s,L)=>s+L.e,0), a=lev.filter(L=>L.k==='a').reduce((s,L)=>s+L.e,0), n=lev.filter(L=>L.k==='n').reduce((s,L)=>s+L.e,0);
  const cnt=(t,k)=>lev.filter(L=>L.t===t&&L.k===k).reduce((s,L)=>s+L.e,0);
  return {lev,ne,bo:(b-a)/2,b,a,n,unp,sb:cnt('σ','b'),pb:cnt('π','b'),sa:cnt('σ','a'),pa:cnt('π','a'),mag:unp?'paramagnético':'diamagnético'}; }
/* ---------- color y fotones ---------- */
const COLW=[[380,435,'violeta','amarillo verdoso'],[435,480,'azul','amarillo'],[480,490,'azul verdoso','naranja'],[490,500,'verde azulado','rojo'],[500,560,'verde','púrpura'],[560,580,'amarillo verdoso','violeta'],[580,595,'amarillo','azul'],[595,605,'naranja','azul verdoso'],[605,750,'rojo','verde azulado']];
function complement(nm){ const r=COLW.find(w=>nm>=w[0]&&nm<w[1]); return r?{abs:r[2],seen:r[3]}:null; }
function photonsToHeat(lamAbs,lamEm,mass_g,dT,cp=4.184){ const dE=h*c*(1/(lamAbs*1e-9)-(lamEm?1/(lamEm*1e-9):0)); const q=mass_g*cp*dT; return {dE,q,n:q/dE,mol:q/dE/NA}; }
return {h,c,NA,eV,RH,rydberg,ionization,region,SERIES,comb,shellStates,termsFromShells,hundGround,parseConfig,J2s,landeG,zeeman,octSpin,operators,LL,
  decayGuess,seriesCount,YEAR,decay,halfLifeFromActivity,bindingEnergy,NUC,spinStates,parseFormula,cluster,fragments,ISO,diatomic,complement,photonsToHeat,COLW,AR,ZOF};
})();
if(typeof module!=='undefined') module.exports=ATOM;
