/* ===== Cálculos avanzados puros (sin DOM) ===== */
const ADV = (function(){
"use strict";
// ---------- espín–órbita: un electrón d en campo octaédrico + ζ l·s ----------
// base |m, σ>, m = 2..−2, σ = +½, −½  (índice 2·im + is)
function soMatrix(Dq, zeta){
  const M=[2,1,0,-1,-2], N=10, H=Array.from({length:N},()=>new Array(N).fill(0));
  const v=(a,b)=>{ if(a===b) return a===0?6:(Math.abs(a)===1?-4:1); if(Math.abs(a)===2&&a===-b) return 5; return 0; };
  const idx=(im,is)=>2*im+is;
  for(let i=0;i<5;i++)for(let j=0;j<5;j++)for(let s=0;s<2;s++) H[idx(i,s)][idx(j,s)]+=Dq*v(M[i],M[j]);
  for(let i=0;i<5;i++){ const m=M[i];
    for(let s=0;s<2;s++){ const ms=s===0?0.5:-0.5; H[idx(i,s)][idx(i,s)]+=zeta*m*ms; }
    // ½ζ(l₊s₋ + l₋s₊): l₊s₋ lleva |m,α> → |m+1,β>
    if(m<2){ const j=M.indexOf(m+1); const c=0.5*zeta*Math.sqrt(6-m*(m+1)); H[idx(j,1)][idx(i,0)]+=c; H[idx(i,0)][idx(j,1)]+=c; }
  }
  return H;
}
function soLevels(Dq,zeta){ return LF.jacobi(soMatrix(Dq,zeta)).vals.slice().sort((a,b)=>a-b); }

// ---------- lantánidos: reglas de Hund ----------
const LN=['La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'];
// Rangos experimentales típicos a 300 K (μ_B). Fuentes: apuntes de magnetismo de SIU (Atkins y Friedman) para Tb–Lu;
// tabla de propiedades de Ln³⁺ (resumen tipo Greenwood) para Ce–Er. Pm: sin dato fiable.
const LNEXP=[[0,0],[2.3,2.5],[3.4,3.6],[3.5,3.6],null,[1.4,1.7],[3.3,3.5],[7.9,8.0],[9.5,9.8],[10.4,10.6],[10.4,10.7],[9.4,9.6],[7.3,7.3],[4.5,4.5],[0,0]];
function hund(n){ // 4f^n
  const ml=[3,2,1,0,-1,-2,-3]; let S=0,L=0;
  for(let k=0;k<n;k++){ if(k<7){S+=0.5;L+=ml[k];} else {S-=0.5;L+=ml[k-7];} }
  const J= n<7 ? Math.abs(L-S) : L+S;
  const g= J===0?0: 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));
  return {n,S,L,J,g,mu:J===0?0:g*Math.sqrt(J*(J+1)),muSO:Math.sqrt(4*S*(S+1))};
}

// ---------- energía de activación de campo cristalino (Dq) ----------
const EOCT=[6,6,-4,-4,-4], ESPY=[0.86,9.14,-0.86,-4.57,-4.57], EPBP=[4.93,2.82,2.82,-5.28,-5.28];
function fill(E,n,P){ let best=null; for(let c=0;c<243;c++){const o=[];let t=c,s=0;for(let i=0;i<5;i++){o.push(t%3);s+=t%3;t=Math.floor(t/3);} if(s!==n)continue; let e=0,p=0,u=0;for(let i=0;i<5;i++){e+=o[i]*E[i];if(o[i]===2)p++;if(o[i]===1)u++;} const tot=e+P*p; if(!best||tot<best.tot-1e-9||(Math.abs(tot-best.tot)<1e-9&&u>best.u))best={tot,e,u,o};} return best; }
function cfae(n,ls){ const P=ls?0:1e9; const o=fill(EOCT,n,P).e, s=fill(ESPY,n,P).e, p=fill(EPBP,n,P).e; return {oct:o,spy:s,pbp:p,dSP:s-o,dPBP:p-o}; }

// ---------- EPR: tensor g de Cu(II) y espectro en polvo ----------
function gCu(o){ const ge=2.0023, l=o.lambda, k2=o.k2;
  if(o.ground==='x2y2') return {gpar:ge-8*l*k2/o.dxy, gperp:ge-2*l*k2/o.dxz};
  return {gpar:ge, gperp:ge-6*l*k2/o.dxz}; }
function eprPowder(o){ // o: gpar,gperp,Apar,Aperp (1e-4 cm-1), nu GHz, I (espín nuclear), lw mT, Bmin,Bmax
  const NB=900, Bmin=o.Bmin, Bmax=o.Bmax, dB=(Bmax-Bmin)/NB, abs=new Float64Array(NB+1);
  const muB=0.46686; // cm-1/T
  const mIs=[]; for(let m=-o.I;m<=o.I+1e-9;m+=1) mIs.push(m);
  const Nt=1500;
  for(let k=0;k<Nt;k++){ const c=(k+0.5)/Nt, s2=1-c*c;
    const g=Math.sqrt(o.gpar*o.gpar*c*c+o.gperp*o.gperp*s2);
    const A=Math.sqrt(o.Apar*o.Apar*o.gpar*o.gpar*c*c+o.Aperp*o.Aperp*o.gperp*o.gperp*s2)/g*1e-4; // cm-1
    const B0=71.4477*o.nu/g; // mT
    const Ab=A/(g*muB)*1000; // mT
    mIs.forEach(m=>{ const B=B0-m*Ab; const i0=Math.round((B-Bmin)/dB);
      for(let i=Math.max(0,i0-40);i<=Math.min(NB,i0+40);i++){ const x=(Bmin+i*dB-B)/o.lw; abs[i]+=Math.exp(-x*x); } });
  }
  const der=new Float64Array(NB+1); for(let i=1;i<NB;i++) der[i]=(abs[i+1]-abs[i-1])/(2*dB);
  return {B:Array.from({length:NB+1},(_,i)=>Bmin+i*dB),abs,der,Bpar:71.4477*o.nu/o.gpar,Bperp:71.4477*o.nu/o.gperp};
}

// ---------- estructuras cristalinas (coordenadas fraccionarias) ----------
function spinel(){ // grupo Fd-3m, origen 2, u ideal = 1/4
  const fcc=[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]], w=x=>((x%1)+1)%1;
  const add=(base,arr,type)=>base.forEach(b=>fcc.forEach(t=>arr.push({p:[w(b[0]+t[0]),w(b[1]+t[1]),w(b[2]+t[2])],t:type})));
  const at=[]; add([[1/8,1/8,1/8],[7/8,3/8,3/8]],at,'A'); add([[.5,.5,.5],[.25,.75,0],[.75,0,.25],[0,.25,.75]],at,'B');
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)for(let k=0;k<4;k++){ const odd=(i%2)+(j%2)+(k%2); if(odd===3||odd===1) at.push({p:[i/4,j/4,k/4],t:'O'}); }
  return at;
}
function perovskite(){ const at=[]; for(const x of [0,1])for(const y of [0,1])for(const z of [0,1]) at.push({p:[x,y,z],t:'A'});
  at.push({p:[.5,.5,.5],t:'B'}); [[.5,.5,0],[.5,.5,1],[.5,0,.5],[.5,1,.5],[0,.5,.5],[1,.5,.5]].forEach(p=>at.push({p,t:'O'})); return at; }
function rocksalt(){ const at=[]; for(let i=0;i<=2;i++)for(let j=0;j<=2;j++)for(let k=0;k<=2;k++){ const p=[i/2,j/2,k/2]; at.push({p,t:((i+j+k)%2===0)?'B':'O'}); } return at; }
function bonds(at,cut,typesA,typesB){ const out=[]; for(let i=0;i<at.length;i++)for(let j=i+1;j<at.length;j++){ const a=at[i],b=at[j];
  if(!((typesA.includes(a.t)&&typesB.includes(b.t))||(typesA.includes(b.t)&&typesB.includes(a.t)))) continue;
  const d=Math.hypot(a.p[0]-b.p[0],a.p[1]-b.p[1],a.p[2]-b.p[2]); if(d<cut) out.push([i,j]); } return out; }
// número de coordinación con condiciones periódicas (para comprobar)
function coordPBC(at,i,type,r){ const a=at[i]; let n=0; at.forEach(b=>{ if(b.t!==type) return; let d2=0; for(let k=0;k<3;k++){ let d=b.p[k]-a.p[k]; d-=Math.round(d); d2+=d*d; } if(Math.abs(Math.sqrt(d2)-r)<1e-6) n++; }); return n; }

// ---------- superficie Jahn–Teller ----------
function jtE(rho,phi,F,K,G){ return 0.5*K*rho*rho - F*rho - G*rho*rho*rho*Math.cos(3*phi); }

// ---------- detección de picos ----------
function parseSpectrum(txt){ const pts=[]; txt.split(/\r?\n/).forEach(l=>{ const m=l.trim().split(/[\s,;\t]+/).map(s=>parseFloat(s.replace(',','.'))); if(m.length>=2&&isFinite(m[0])&&isFinite(m[1])) pts.push([m[0],m[1]]); });
  if(!pts.length) return null; const nm=pts.every(p=>p[0]<2500); const cm=pts.map(p=>[nm?1e7/p[0]:p[0],p[1]]).sort((a,b)=>a[0]-b[0]); return {pts:cm,unit:nm?'nm':'cm⁻¹'}; }
function findPeaks(pts,minRel=0.06){ const y=pts.map(p=>p[1]); const sm=y.map((v,i)=>{let s=0,c=0;for(let k=-3;k<=3;k++){if(y[i+k]!=null){s+=y[i+k];c++;}}return s/c;});
  const mx=Math.max(...sm), mn=Math.min(...sm), pk=[];
  for(let i=2;i<sm.length-2;i++){ if(sm[i]>=sm[i-1]&&sm[i]>=sm[i+1]&&sm[i]>sm[i-2]&&sm[i]>sm[i+2]&&sm[i]-mn>minRel*(mx-mn)){
      // prominencia: descenso mínimo a ambos lados
      let l=sm[i],r=sm[i]; for(let k=i;k>=0;k--){l=Math.min(l,sm[k]); if(sm[k]>sm[i])break;} for(let k=i;k<sm.length;k++){r=Math.min(r,sm[k]); if(sm[k]>sm[i])break;}
      if(sm[i]-Math.max(l,r)>minRel*(mx-mn)) pk.push({x:pts[i][0],y:sm[i]}); } }
  return pk; }
return {soLevels,soMatrix,hund,LN,LNEXP,cfae,gCu,eprPowder,spinel,perovskite,rocksalt,bonds,coordPBC,jtE,parseSpectrum,findPeaks};
})();
if(typeof module!=='undefined') module.exports=ADV;
