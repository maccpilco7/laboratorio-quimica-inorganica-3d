m={'ENGINE':'engine.js','ADV':'adv_pure.js','GLOS':'glos.js','CRY':'cryst.js','MOL':'mol.js','ATOM':'atom.js','EXAMS':'exams.js','PRAC':'prac.js','INORG':'inorg.js','APP':'app_full.js','SYM':'sym.js','SHADERS':'shaders.js','LAB3D':'lab3d.js'}
h=open('src.html').read()
for k,f in m.items():
    h=h.replace('/*'+k+'*/',open(f).read().replace('module.exports','module.exports'),1)
open('index.html','w').write(h); print(len(h))
