/* ===== Shaders GLSL ES 3.00 para raymarching de isosuperficies =====
   Se usan con THREE.RawShaderMaterial + glslVersion: THREE.GLSL3 (Three antepone "#version 300 es").
   Campo: ψ(p) = Σ cᵢ φᵢ(M⁻¹ p) — la operación de simetría se aplica al argumento (representación activa). */
const ORB_SHADERS = {
vert: `
precision highp float;
uniform mat4 modelMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
in vec3 position;
out vec3 vWorld;
void main(){
  vWorld = (modelMatrix * vec4(position, 1.0)).xyz;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`,
frag: `
precision highp float;
precision highp int;
uniform vec3  uCam;          // posición de la cámara (mundo = coordenadas químicas, z arriba)
uniform mat4  uViewProj;     // para escribir gl_FragDepth y ocluir correctamente con los ligandos
uniform mat3  uMinv;         // inversa de la matriz de la operación (interpolada) × orientación acumulada
uniform float uC[17];        // coeficientes: 5 d, 3 p, s, 8 ligandos (SALC)
uniform vec3  uLig[8];       // posiciones de los ligandos (base, sin transformar)
uniform vec3  uPiV[8];       // coeficiente vectorial del orbital p de cada ligando (SALCs π)
uniform float uIso;          // isovalor |ψ|
uniform float uOpacity;
uniform vec3  uColPos;
uniform vec3  uColNeg;
uniform float uBox;          // semiarista de la caja de marcha
uniform int   uSteps;
in vec3 vWorld;
out vec4 fragColor;

// ψ evaluada en coordenadas de la función base
float field(vec3 p){
  vec3 q = uMinv * p;
  float r2 = dot(q, q);
  float r  = sqrt(r2);
  // parte radial hidrogenoide 3d: ρ² e^(−ρ), ρ = 1.6 r ; normalizada para que max|d_z²| = 1
  float rad = exp(-1.6 * r) * 2.56 / (0.5413 * 1.118);
  float d = uC[0] * 0.5590 * (3.0*q.z*q.z - r2)
          + uC[1] * 0.9682 * (q.x*q.x - q.y*q.y)
          + uC[2] * 1.9365 * q.x*q.y
          + uC[3] * 1.9365 * q.x*q.z
          + uC[4] * 1.9365 * q.y*q.z;
  float v = d * rad;
  float pr = 1.6 * exp(-1.6 * r) / 0.3679;              // 4p: ρ e^(−ρ)
  v += (uC[5]*q.x + uC[6]*q.y + uC[7]*q.z) * pr;
  v += uC[8] * exp(-r2 / 0.6);                          // 4s (gaussiana)
  for(int k = 0; k < 8; k++){                           // orbitales σ de los ligandos
    vec3 dl = q - uLig[k];
    float gk = exp(-dot(dl, dl) / 0.32);
    v += 1.5 * uC[9 + k] * gk;                          // σ: función tipo s
    v += 6.18 * dot(dl, uPiV[k]) * gk;                  // π: función tipo p (1,5 × 1/0,2426)
  }
  return v;
}
float sfield(vec3 p){ return abs(field(p)) - uIso; }
vec3 gradAbs(vec3 p){
  const float e = 0.012;
  return vec3(abs(field(p+vec3(e,0,0))) - abs(field(p-vec3(e,0,0))),
              abs(field(p+vec3(0,e,0))) - abs(field(p-vec3(0,e,0))),
              abs(field(p+vec3(0,0,e))) - abs(field(p-vec3(0,0,e))));
}
vec2 boxHit(vec3 ro, vec3 rd){
  vec3 inv = 1.0 / rd;
  vec3 t0 = (vec3(-uBox) - ro) * inv, t1 = (vec3(uBox) - ro) * inv;
  vec3 tmin = min(t0, t1), tmax = max(t0, t1);
  return vec2(max(max(tmin.x, tmin.y), tmin.z), min(min(tmax.x, tmax.y), tmax.z));
}
vec3 shade(vec3 p, vec3 rd, float phase){
  vec3 n = -normalize(gradAbs(p) + 1e-9);               // normal saliente de la isosuperficie
  if(dot(n, rd) > 0.0) n = -n;                          // dos caras
  vec3 base = phase > 0.0 ? uColPos : uColNeg;
  vec3 L1 = normalize(vec3(0.45, -0.35, 0.82));
  vec3 V = -rd;
  float dif = 0.62*max(dot(n, L1), 0.0) + 0.38*max(dot(n, V), 0.0);
  vec3 H = normalize(L1 + V);
  float spec = pow(max(dot(n, H), 0.0), 48.0);
  float fres = pow(1.0 - max(dot(n, V), 0.0), 3.0);     // aproximación de Schlick
  return base*(0.16 + dif) + vec3(0.55)*spec + base*fres*0.45;
}
void main(){
  vec3 ro = uCam;
  vec3 rd = normalize(vWorld - uCam);
  vec2 th = boxHit(ro, rd);
  float t0 = max(th.x, 0.0), t1 = th.y;
  if(t1 <= t0) discard;
  float dt = (t1 - t0) / float(uSteps);
  float t = t0;
  float prev = sfield(ro + rd*t);
  vec4 acc = vec4(0.0);
  float tFirst = -1.0;
  for(int i = 0; i < 400; i++){
    if(i >= uSteps) break;
    t += dt;
    float cur = sfield(ro + rd*t);
    if((prev < 0.0) != (cur < 0.0)){
      float a = t - dt, b = t;
      for(int j = 0; j < 6; j++){                       // refinamiento por bisección
        float m = 0.5*(a + b);
        if((sfield(ro + rd*m) < 0.0) == (prev < 0.0)) a = m; else b = m;
      }
      vec3 h = ro + rd*b;
      vec3 c = shade(h, rd, sign(field(h)));
      acc.rgb += (1.0 - acc.a) * uOpacity * c;
      acc.a   += (1.0 - acc.a) * uOpacity;
      if(tFirst < 0.0) tFirst = b;
      if(acc.a > 0.97) break;
    }
    prev = cur;
  }
  if(acc.a < 0.01) discard;
  fragColor = vec4(acc.rgb / acc.a, acc.a);
  vec4 clip = uViewProj * vec4(ro + rd*tFirst, 1.0);
  gl_FragDepth = clamp(clip.z/clip.w*0.5 + 0.5, 0.0, 1.0);
}`
};
if(typeof module!=='undefined') module.exports=ORB_SHADERS;
