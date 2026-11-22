/* ===== Motor numérico: CI d^n (Tanabe–Sugano), AOM, autovalores ===== */
const LF = (function(){
"use strict";
// ---------- álgebra ----------
function jacobi(Ain){ // simétrica pequeña, devuelve {vals, vecs(cols)}
  const n=Ain.length, A=Ain.map(r=>r.slice()), V=[];
  for(let i=0;i<n;i++){V.push(new Array(n).fill(0));V[i][i]=1;}
  for(let sweep=0;sweep<60;sweep++){
    let off=0; for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)off+=A[i][j]*A[i][j];
    if(off<1e-22) break;
    for(let p=0;p<n;p++)for(let q=p+1;q<n;q++){
      if(Math.abs(A[p][q])<1e-14) continue;
      const th=(A[q][q]-A[p][p])/(2*A[p][q]);
      const t=Math.sign(th||1)/(Math.abs(th)+Math.sqrt(th*th+1));
      const c=1/Math.sqrt(t*t+1), s=t*c;
      for(let k=0;k<n;k++){const akp=A[k][p],akq=A[k][q];A[k][p]=c*akp-s*akq;A[k][q]=s*akp+c*akq;}
      for(let k=0;k<n;k++){const apk=A[p][k],aqk=A[q][k];A[p][k]=c*apk-s*aqk;A[q][k]=s*apk+c*aqk;}
      for(let k=0;k<n;k++){const vkp=V[k][p],vkq=V[k][q];V[k][p]=c*vkp-s*vkq;V[k][q]=s*vkp+c*vkq;}
    }
  }
  return {vals:A.map((r,i)=>r[i]), vecs:V};
}
// Householder + QL implícito (solo autovalores)
function eigvals(Ain){
  const n=Ain.length; if(n===1) return [Ain[0][0]];
  const a=Ain.map(r=>Float64Array.from(r)), d=new Float64Array(n), e=new Float64Array(n);
  for(let i=n-1;i>0;i--){
    const l=i-1; let h=0,scale=0;
    if(l>0){
      for(let k=0;k<=l;k++) scale+=Math.abs(a[i][k]);
      if(scale===0) e[i]=a[i][l];
      else{
        for(let k=0;k<=l;k++){a[i][k]/=scale;h+=a[i][k]*a[i][k];}
        let f=a[i][l]; let g=f>=0?-Math.sqrt(h):Math.sqrt(h);
        e[i]=scale*g; h-=f*g; a[i][l]=f-g; f=0;
        for(let j=0;j<=l;j++){
          g=0; for(let k=0;k<=j;k++) g+=a[j][k]*a[i][k];
          for(let k=j+1;k<=l;k++) g+=a[k][j]*a[i][k];
          e[j]=g/h; f+=e[j]*a[i][j];
        }
        const hh=f/(h+h);
        for(let j=0;j<=l;j++){
          f=a[i][j]; e[j]=g=e[j]-hh*f;
          for(let k=0;k<=j;k++) a[j][k]-=(f*e[k]+g*a[i][k]);
        }
      }
    } else e[i]=a[i][l];
    d[i]=h;
  }
  for(let i=0;i<n;i++) d[i]=a[i][i];
  for(let i=1;i<n;i++) e[i-1]=e[i]; e[n-1]=0;
  for(let l=0;l<n;l++){
    let iter=0,m;
    do{
      for(m=l;m<n-1;m++){const dd=Math.abs(d[m])+Math.abs(d[m+1]); if(Math.abs(e[m])<=1e-15*dd) break;}
      if(m!==l){
        if(iter++>60) break;
        let g=(d[l+1]-d[l])/(2*e[l]); let r=Math.hypot(g,1);
        g=d[m]-d[l]+e[l]/(g+(g>=0?Math.abs(r):-Math.abs(r)));
        let s=1,c=1,p=0,i;
        for(i=m-1;i>=l;i--){
          let f=s*e[i]; const b=c*e[i];
          e[i+1]=(r=Math.hypot(f,g));
          if(r===0){d[i+1]-=p;e[m]=0;break;}
          s=f/r;c=g/r;g=d[i+1]-p;r=(d[i]-g)*s+2*c*b;d[i+1]=g+(p=s*r);g=c*r-b;
        }
        if(r===0&&i>=l) continue;
        d[l]-=p;e[l]=g;e[m]=0;
      }
    }while(m!==l);
  }
  return Array.from(d).sort((x,y)=>x-y);
}

// ---------- coeficientes c^k (Condon–Shortley) ----------
const fact=[1]; for(let i=1;i<30;i++) fact[i]=fact[i-1]*i;
function w3j(j1,j2,j3,m1,m2,m3){
  if(m1+m2+m3!==0) return 0;
  if(j3<Math.abs(j1-j2)||j3>j1+j2) return 0;
  if(Math.abs(m1)>j1||Math.abs(m2)>j2||Math.abs(m3)>j3) return 0;
  const tri=fact[j1+j2-j3]*fact[j1-j2+j3]*fact[-j1+j2+j3]/fact[j1+j2+j3+1];
  const pre=Math.sqrt(tri*fact[j1+m1]*fact[j1-m1]*fact[j2+m2]*fact[j2-m2]*fact[j3+m3]*fact[j3-m3]);
  let s=0;
  for(let t=0;t<=40;t++){
    const a=[t, j3-j2+t+m1, j3-j1+t-m2, j1+j2-j3-t, j1-t-m1, j2-t+m2];
    if(a.some(x=>x<0)) continue;
    s+=(t%2?-1:1)/(fact[a[0]]*fact[a[1]]*fact[a[2]]*fact[a[3]]*fact[a[4]]*fact[a[5]]);
  }
  return ((j1-j2-m3)%2?-1:1)*pre*s;
}
const M=[2,1,0,-1,-2];
function ck(k,m,mp){ return ((m%2)?-1:1)*5*w3j(2,k,2,0,0,0)*w3j(2,k,2,-m,m-mp,mp); }
// <ab|cd> para un k dado
function I2(k,a,b,c,d){ if(a+b!==c+d) return 0; return ck(k,a,c)*ck(k,d,b); }
// potencial octaédrico en base |m> (unidades Dq)
function vOct(mi,mj){ const a=M[mi],b=M[mj];
  if(a===b) return a===0?6:(Math.abs(a)===1?-4:1);
  if(Math.abs(a)===2&&a===-b) return 5; return 0; }

// ---------- bloques de determinantes ----------
const popc=x=>{let c=0;while(x){c+=x&1;x>>=1;}return c;};
function combos(n,k){const out=[];(function rec(s,start,cur){if(cur.length===k){out.push(cur.slice());return;}for(let i=start;i<n;i++){cur.push(i);rec(s,i+1,cur);cur.pop();}})(0,0,[]);return out;}
const blockCache={};
function block(n,nu){ // nu electrones α, n-nu β
  const key=n+'_'+nu; if(blockCache[key]) return blockCache[key];
  const nd=n-nu; const dets=[];
  combos(5,nu).forEach(U=>combos(5,nd).forEach(Dn=>{let m=0;U.forEach(i=>m|=1<<(2*i));Dn.forEach(i=>m|=1<<(2*i+1));dets.push(m);}));
  const idx=new Map(dets.map((d,i)=>[d,i])); const N=dets.length;
  const mk=()=>Array.from({length:N},()=>new Float64Array(N));
  const V=mk(),G2=mk(),G4=mk();
  const sgn=(mask,p)=>(popc(mask&((1<<p)-1))%2)?-1:1;
  dets.forEach((D,j)=>{
    // un electrón
    for(let p=0;p<10;p++){ if(!(D>>p&1)) continue;
      for(let r=p%2;r<10;r+=2){ const v=vOct(r>>1,p>>1); if(!v) continue;
        if(r===p){V[j][j]+=v;continue;}
        if(D>>r&1) continue;
        let s=sgn(D,p); let D1=D&~(1<<p); s*=sgn(D1,r); const D2=D1|(1<<r);
        V[idx.get(D2)][j]+=s*v; } }
    // dos electrones
    for(let p=0;p<10;p++){ if(!(D>>p&1)) continue;
      for(let q=p+1;q<10;q++){ if(!(D>>q&1)) continue;
        let s1=sgn(D,p); let D1=D&~(1<<p); s1*=sgn(D1,q); D1&=~(1<<q);
        for(let r=0;r<10;r++){ if(D1>>r&1) continue;
          for(let s=r+1;s<10;s++){ if(D1>>s&1) continue;
            if((r%2)+(s%2)!==(p%2)+(q%2)) continue;
            const mr=M[r>>1],ms=M[s>>1],mp=M[p>>1],mq=M[q>>1];
            const dir=(r%2===p%2&&s%2===q%2), exc=(r%2===q%2&&s%2===p%2);
            let a2=0,a4=0;
            if(dir){a2+=I2(2,mr,ms,mp,mq);a4+=I2(4,mr,ms,mp,mq);}
            if(exc){a2-=I2(2,mr,ms,mq,mp);a4-=I2(4,mr,ms,mq,mp);}
            if(Math.abs(a2)<1e-14&&Math.abs(a4)<1e-14) continue;
            let sg=s1*sgn(D1,s); const D2=D1|(1<<s); sg*=sgn(D2,r); const D3=D2|(1<<r);
            const i=idx.get(D3); G2[i][j]+=sg*a2; G4[i][j]+=sg*a4;
          } } } }
  });
  return blockCache[key]={N,V,G2,G4};
}
function Hmat(b,dq,B,C){ const F2=49*B+7*C, F4=12.6*C; const N=b.N, H=[];
  for(let i=0;i<N;i++){const row=new Array(N);for(let j=0;j<N;j++)row[j]=dq*b.V[i][j]+F2*b.G2[i][j]+F4*b.G4[i][j];H.push(row);} return H; }
function smax(n){return Math.min(n,10-n)/2;}
// estados por espín: {S: [energías ordenadas]}
function spectrum(n,Dq,B,C,onlyS){
  const S0=smax(n), lists={};
  for(let S=S0;S>=0;S-=1){ if(onlyS!==undefined && S<onlyS) break; const nu=Math.round(n/2+S); lists[S]=eigvals(Hmat(block(n,nu),Dq,B,C)); }
  const out={};
  for(let S=S0;S>=0;S-=1){ if(!lists[S]) break; const a=lists[S], b=lists[S+1]||[]; const res=[]; let j=0;
    const tol=1e-6*Math.max(1,Math.abs(B)*50);
    for(let i=0;i<a.length;i++){ if(j<b.length && Math.abs(a[i]-b[j])<tol){j++;continue;} res.push(a[i]); }
    out[S]=res; }
  return out;
}
function groupLevels(arr,tol){ const g=[]; arr.forEach(e=>{const last=g[g.length-1]; if(last&&Math.abs(last.e-e)<tol){last.n++;} else g.push({e,n:1});}); return g; }


// ---------- ajuste de bandas (estados de espín máximo) ----------
function transitions(n,Dq,B,C){ const S=smax(n); const lv=groupLevels(spectrum(n,Dq,B,C,S)[S],1); return lv.slice(1).map(l=>({e:l.e-lv[0].e,n:l.n})); }
function bestAssign(obs,calc){ let best=Infinity, pick=null; const k=obs.length, m=calc.length;
  (function rec(i,start,acc,sel){ if(i===k){ if(acc<best){best=acc;pick=sel.slice();} return; }
    for(let j=start;j<=m-(k-i);j++){ const d=obs[i]-calc[j].e; sel.push(j); rec(i+1,j+1,acc+d*d,sel); sel.pop(); } })(0,0,0,[]);
  return {sse:best,pick}; }
function fit(n,obs,gamma){ obs=obs.slice().sort((a,b)=>a-b); let best={sse:Infinity};
  const evalAt=(D,B)=>{ const c=transitions(n,D/10,B,gamma*B); if(c.length<obs.length) return; const r=bestAssign(obs,c); if(r.sse<best.sse) best={sse:r.sse,D,B,calc:c,pick:r.pick}; };
  for(let D=5000;D<=36000;D+=500) for(let B=300;B<=1200;B+=50) evalAt(D,B);
  const D0=best.D,B0=best.B;
  for(let D=D0-600;D<=D0+600;D+=20) for(let B=Math.max(200,B0-60);B<=B0+60;B+=4) evalAt(D,B);
  best.rms=Math.sqrt(best.sse/obs.length); best.obs=obs; return best; }

// ---------- AOM ----------
const r3=Math.sqrt(3)/2;
const AM=[ [[-.5,0,0],[0,-.5,0],[0,0,1]], [[r3,0,0],[0,-r3,0],[0,0,0]], [[0,r3,0],[r3,0,0],[0,0,0]], [[0,0,r3],[0,0,0],[r3,0,0]], [[0,0,0],[0,0,r3],[0,r3,0]] ];
const qf=(A,u,v)=>{let s=0;for(let i=0;i<3;i++)for(let j=0;j<3;j++)s+=u[i]*A[i][j]*v[j];return s;};
function perp(n){ const a=Math.abs(n[2])<0.9?[0,0,1]:[1,0,0];
  let u=[a[1]*n[2]-a[2]*n[1],a[2]*n[0]-a[0]*n[2],a[0]*n[1]-a[1]*n[0]]; const l=Math.hypot(...u); u=u.map(x=>x/l);
  const v=[n[1]*u[2]-n[2]*u[1],n[2]*u[0]-n[0]*u[2],n[0]*u[1]-n[1]*u[0]]; return [u,v]; }
function norm(v){const l=Math.hypot(...v);return v.map(x=>x/l);}
function aomMatrix(sites){ // sites: [{n:[x,y,z], es, ep}]
  const H=[0,1,2,3,4].map(()=>[0,0,0,0,0]);
  sites.forEach(s=>{ const n=norm(s.n); const [u,v]=perp(n);
    const Fs=AM.map(A=>qf(A,n,n)), Fu=AM.map(A=>qf(A,n,u)/r3), Fv=AM.map(A=>qf(A,n,v)/r3);
    for(let i=0;i<5;i++)for(let j=0;j<5;j++) H[i][j]+=s.es*Fs[i]*Fs[j]+s.ep*(Fu[i]*Fu[j]+Fv[i]*Fv[j]); });
  return H;
}
function aomSolve(sites){ const r=jacobi(aomMatrix(sites)); const ord=[0,1,2,3,4].sort((a,b)=>r.vals[a]-r.vals[b]);
  return ord.map(k=>({e:r.vals[k], v:[0,1,2,3,4].map(i=>r.vecs[i][k])})); }
const deg=Math.PI/180;
function sph(th,ph){return [Math.sin(th*deg)*Math.cos(ph*deg),Math.sin(th*deg)*Math.sin(ph*deg),Math.cos(th*deg)];}
const GEOM={
  oct:{name:'Octaédrica',pts:[[0,0,1],[0,0,-1],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]]},
  tet:{name:'Tetraédrica',pts:[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(norm)},
  sqp:{name:'Plano cuadrada',pts:[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]]},
  spy:{name:'Pirámide cuadrada',pts:[[0,0,1],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]]},
  tbp:{name:'Bipirámide trigonal',pts:[[0,0,1],[0,0,-1],sph(90,0),sph(90,120),sph(90,240)]},
  tpr:{name:'Prisma trigonal',pts:[sph(49.1,0),sph(49.1,120),sph(49.1,240),sph(130.9,0),sph(130.9,120),sph(130.9,240)]},
  pbp:{name:'Bipirámide pentagonal',pts:[[0,0,1],[0,0,-1],sph(90,0),sph(90,72),sph(90,144),sph(90,216),sph(90,288)]},
  cub:{name:'Cúbica',pts:[[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]].map(norm)},
  trp:{name:'Trigonal plana',pts:[sph(90,0),sph(90,120),sph(90,240)]},
  lin:{name:'Lineal',pts:[[0,0,1],[0,0,-1]]}
};
return {fit,transitions,jacobi,eigvals,spectrum,groupLevels,smax,block,Hmat,aomSolve,aomMatrix,GEOM,sph,norm};
})();
if(typeof module!=='undefined') module.exports=LF;
