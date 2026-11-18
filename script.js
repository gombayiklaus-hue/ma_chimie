/* ---------- Données (éléments 1 à 54) ---------- */
const SYMS="H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe".split(",");
const NAMES="Hydrogène,Hélium,Lithium,Béryllium,Bore,Carbone,Azote,Oxygène,Fluor,Néon,Sodium,Magnésium,Aluminium,Silicium,Phosphore,Soufre,Chlore,Argon,Potassium,Calcium,Scandium,Titane,Vanadium,Chrome,Manganèse,Fer,Cobalt,Nickel,Cuivre,Zinc,Gallium,Germanium,Arsenic,Sélénium,Brome,Krypton,Rubidium,Strontium,Yttrium,Zirconium,Niobium,Molybdène,Technétium,Ruthénium,Rhodium,Palladium,Argent,Cadmium,Indium,Étain,Antimoine,Tellure,Iode,Xénon".split(",");
// Exceptions à la règle de Klechkowski
const EXC={24:{"4s":1,"3d":5},29:{"4s":1,"3d":10},41:{"5s":1,"4d":4},42:{"5s":1,"4d":5},44:{"5s":1,"4d":7},45:{"5s":1,"4d":8},46:{"5s":0,"4d":10},47:{"5s":1,"4d":10}};
const NG=[[2,"He"],[10,"Ne"],[18,"Ar"],[36,"Kr"],[54,"Xe"]];
// Géométries VSEPR : clé "X-E" -> [nom, angle]
const GEO={"2-0":["Linéaire","180°"],"3-0":["Triangulaire plane","120°"],"2-1":["Coudée (angulaire)","≈ 119°"],"4-0":["Tétraédrique","109,5°"],"3-1":["Pyramide trigonale","≈ 107°"],"2-2":["Coudée (angulaire)","≈ 104,5°"],"5-0":["Bipyramide trigonale","90° et 120°"],"4-1":["Bascule (tétraèdre déformé)","≈ 90° et 120°"],"3-2":["En T","≈ 90°"],"2-3":["Linéaire","180°"],"6-0":["Octaédrique","90°"],"5-1":["Pyramide à base carrée","≈ 90°"],"4-2":["Carré plan","90°"],"7-0":["Bipyramide à base pentagonale","72° et 90°"]};
const BO={H:1,F:1,Cl:1,Br:1,I:1,O:2,S:2,Se:2,N:3}; // liaisons d'un atome terminal

/* ---------- Structure électronique ---------- */
const $=id=>document.getElementById(id);
const key=(n,l)=>n+"spdf"[l], cap=l=>4*l+2, L=k=>"spdf".indexOf(k[1]);
const ORDER=[];for(let n=1;n<=6;n++)for(let l=0;l<Math.min(n,4);l++)ORDER.push([n,l]);
ORDER.sort((a,b)=>(a[0]+a[1])-(b[0]+b[1])||a[0]-b[0]); // règle de Klechkowski
function fill(N){const o={};let r=N;for(const[n,l]of ORDER){if(r<=0)break;const e=Math.min(cap(l),r);o[key(n,l)]=e;r-=e;}return o;}
const neutral=Z=>Object.assign(fill(Z),EXC[Z]||{});
function ion(Z,q){const o=neutral(Z);
  for(let i=0;i<q;i++){const ks=Object.keys(o).filter(k=>o[k]>0).sort((a,b)=>b[0]-a[0]||L(b)-L(a));o[ks[0]]--;} // cation : on retire d'abord la couche n la plus haute
  for(let i=0;i<-q;i++)for(const[n,l]of ORDER){const k=key(n,l);if((o[k]||0)<cap(l)){o[k]=(o[k]||0)+1;break;}} // anion : on remplit
  return o;}
const sup=n=>String(n).replace(/\d/g,d=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[d]);
const subStr=s=>String(s).replace(/\d/g,d=>"₀₁₂₃₄₅₆₇₈₉"[d]);
const sub=n=>n>1?subStr(n):"";
const fmt=o=>Object.keys(o).filter(k=>o[k]>0).sort((a,b)=>a[0]-b[0]||L(a)-L(b)).map(k=>k+sup(o[k])).join(" ");
function core(o){for(const[z,s]of NG.slice().reverse()){const f=fill(z);
  if(Object.keys(f).every(k=>o[k]===f[k])){const r={};for(const k in o)if(!(k in f)&&o[k]>0)r[k]=o[k];if(Object.keys(r).length)return`[${s}] ${fmt(r)}`;}}return null;}
function position(Z){const f=fill(Z),o=neutral(Z);
  const per=Math.max(...Object.keys(f).map(k=>+k[0]));
  const occ=ORDER.filter(([n,l])=>(o[key(n,l)]||0)>0),bl="spdf"[occ[occ.length-1][1]];
  const s=o[key(per,0)]||0;let g,ve;
  if(bl==="s"){g=Z===2?18:s;ve=s;}
  else if(bl==="p"){ve=s+(o[key(per,1)]||0);g=10+ve;}
  else{g=s+(o[key(per-1,2)]||0);ve=s;}
  return{per,bl,g,ve};}
const V=s=>position(SYMS.indexOf(s)+1).ve; // électrons de valence
function excitedValence(o){const nmax=Math.max(...Object.keys(o).filter(k=>o[k]>0).map(k=>+k[0]));
  const s=key(nmax,0),p=key(nmax,1);
  if(nmax>1&&o[s]===2&&(o[p]||0)<6){const c={...o};c[s]--;c[p]=(c[p]||0)+1;return c;}
  return null;}
function excitedFirst(o){const occ=ORDER.filter(([n,l])=>(o[key(n,l)]||0)>0),last=occ[occ.length-1];
  for(let j=ORDER.indexOf(last)+1;j<ORDER.length;j++){const[n,l]=ORDER[j],k=key(n,l);
    if((o[k]||0)<cap(l)){const c={...o};c[key(...last)]--;c[k]=(c[k]||0)+1;return c;}}
  return null;}

/* ---------- Dessin de Lewis (SVG) ---------- */
const rad=a=>a*Math.PI/180;
const wrap=s=>`<svg viewBox="0 0 300 220" class="lewis" role="img" aria-label="Schéma de Lewis">${s}</svg>`;
const dot=(cx,cy,r,a,d)=>`<circle cx="${(cx+r*Math.cos(rad(a))-d*Math.sin(rad(a))).toFixed(1)}" cy="${(cy+r*Math.sin(rad(a))+d*Math.cos(rad(a))).toFixed(1)}" r="2.6"/>`;
const pair=(cx,cy,r,a)=>dot(cx,cy,r,a,-5)+dot(cx,cy,r,a,5);
function lonePairs(cx,cy,r,bonds,n){const cand=[0,45,90,135,180,225,270,315],used=[];let o="";
  for(let i=0;i<n;i++){let best=0,bs=-1;
    for(const a of cand){if(used.includes(a))continue;
      const d=Math.min(360,...bonds.concat(used).map(b=>{const x=Math.abs(a-b)%360;return Math.min(x,360-x);}));
      if(d>bs){bs=d;best=a;}}
    used.push(best);o+=pair(cx,cy,r,best);}
  return o;}
function atomLewis(sym,n){let o=`<text x="150" y="110">${sym}</text>`;
  for(let i=0;i<n;i++)o+=dot(150,110,26,[270,0,90,180][i%4],n<=4?0:(i<4?-5:5));
  return wrap(o);}
function drawMol(m){const cx=150,cy=110,R=80,n=m.terms.length;let o="";
  const angs=m.terms.map((_,i)=>n===1?0:n===2?[180,0][i]:(270+i*360/n)%360);
  m.terms.forEach((t,i)=>{const a=rad(angs[i]),px=cx+R*Math.cos(a),py=cy+R*Math.sin(a),k=m.bo[i];
    for(const d of {1:[0],2:[-3.5,3.5],3:[-6,0,6]}[k]){const ox=-d*Math.sin(a),oy=d*Math.cos(a);
      o+=`<line x1="${(cx+18*Math.cos(a)+ox).toFixed(1)}" y1="${(cy+18*Math.sin(a)+oy).toFixed(1)}" x2="${(px-18*Math.cos(a)+ox).toFixed(1)}" y2="${(py-18*Math.sin(a)+oy).toFixed(1)}"/>`;}
    o+=`<text x="${px.toFixed(1)}" y="${py.toFixed(1)}">${t}</text>`+lonePairs(px,py,24,[angs[i]+180],Math.floor((V(t)-k)/2));});
  return wrap(o+`<text x="${cx}" y="${cy}">${m.c}</text>`+lonePairs(cx,cy,24,angs,m.lp));}

/* ---------- Analyse d'un atome ---------- */
const norm=s=>s.trim().normalize("NFD").replace(/[\u0300-\u036f]/g,"").toLowerCase();
function findZ(q){q=norm(q);if(/^\d+$/.test(q)){const z=+q;return z>=1&&z<=54?z:0;}
  let i=SYMS.findIndex(s=>norm(s)===q);if(i<0)i=NAMES.findIndex(n=>norm(n)===q);return i+1;}
const row=(k,v)=>`<div class="row"><dt>${k}</dt><dd>${v}</dd></div>`;
const cfg=o=>{const c=core(o);return`<span class="cfg">${fmt(o)}</span>`+(c?`<br><span class="sub">${c}</span>`:"");};

function showAtom(){const out=$("outAtome"),Z=findZ($("inAtome").value),q=+$("inCharge").value||0;
  if(!Z){out.innerHTML=`<p class="bad">Élément introuvable. Essayez un symbole (Fe), un nom (Carbone) ou un numéro (de 1 à 54).</p>`;return;}
  if(q>=Z){out.innerHTML=`<p class="bad">Charge impossible : ${NAMES[Z-1]} n'a que ${Z} électrons.</p>`;return;}
  const sym=SYMS[Z-1],p=position(Z),base=q?ion(Z,q):neutral(Z);
  const ionName=sym+(q?(Math.abs(q)>1?sup(Math.abs(q)):"")+(q>0?"⁺":"⁻"):"");
  let h=row("Élément",`${NAMES[Z-1]} (${sym})`)+row("Numéro atomique",`Z = ${Z} (${Z} protons)`)
   +row("Position",`Période ${p.per} · groupe ${p.g} · bloc ${p.bl}`)
   +row("État fondamental",cfg(neutral(Z)));
  if(q)h+=row(`Ion ${ionName}`,`${Z-q} électrons<br>${cfg(base)}`);
  const ev=excitedValence(base),ef=excitedFirst(base);
  if(ev)h+=row("Excité (valence)",cfg(ev));
  if(ef)h+=row("1er état excité",cfg(ef));
  h+=row("Électrons de valence",p.bl==="d"?`${p.g} (ns + (n−1)d)`:p.ve)
   +row("Schéma de Lewis",q?"Donné pour l'atome neutre :<br>":"")+row("",atomLewis(sym,p.bl==="d"?p.ve:p.ve))
   +row("Existence",`<span class="ok">✓ Élément connu${q?" ; l'ion "+ionName+" est décrit par la structure ci-dessus":""}</span>`);
  out.innerHTML=h;}

/* ---------- Analyse d'une molécule ---------- */
function analyse(raw){
  const f=raw.replace(/[₀-₉]/g,c=>c.charCodeAt(0)-8320).replace(/\s/g,"");
  const parts=[...f.matchAll(/([A-Z][a-z]?)(\d*)/g)];
  if(!f||parts.map(p=>p[0]).join("")!==f)return{err:"Formule non reconnue. Exemples : H2O, CO2, NH3, SF6."};
  const el={};
  for(const[,s,n]of parts){const z=SYMS.indexOf(s)+1;
    if(!z)return{err:`Élément « ${s} » inconnu ou non géré (Z ≤ 54).`};
    if(position(z).bl==="d")return{err:`${s} est un métal de transition : non géré dans cette version.`};
    el[s]=(el[s]||0)+(+n||1);}
  const cands=Object.keys(el).filter(s=>el[s]===1&&s!=="H").sort((a,b)=>V(a)-V(b)||SYMS.indexOf(b)-SYMS.indexOf(a));
  if(!cands.length)return{err:"Aucun atome central unique : cette version gère les molécules de type AXₙEₘ (un seul atome central)."};
  const c=cands[0],terms=[];
  for(const s in el)for(let i=0;i<(s===c?0:el[s]);i++)terms.push(s);
  if(!terms.length)return{err:"Entrez une molécule (au moins deux atomes)."};
  if(terms.some(t=>!(t in BO)))return{err:"Atome terminal non géré (H, F, Cl, Br, I, O, S, Se, N)."};
  const vc=V(c),per=position(SYMS.indexOf(c)+1).per;
  const total=Object.keys(el).reduce((a,s)=>a+V(s)*el[s],0);
  const bo=terms.map(t=>BO[t]),sum=()=>bo.reduce((a,b)=>a+b,0);
  while(sum()>vc){const i=bo.findIndex(b=>b>1);if(i<0)break;bo[i]--;}
  const base={formula:f,total,c,terms,bo,per,vc,left:vc-sum()};
  if(base.left<0)return{...base,exists:false,why:`${c} n'a que ${vc} électrons de valence : il ne peut pas former ${terms.length} liaisons avec ${terms.join(", ")}.`};
  const lp=Math.floor(base.left/2),rad1=base.left%2,X=terms.length,st=X+lp;
  if(per===2&&st>4)return{...base,lp,exists:false,why:`${c} appartient à la période 2 : pas d'orbitale d, donc 4 doublets au maximum autour de lui (règle de l'octet). Il en faudrait ${st}.`};
  const g=X===1?["Diatomique : linéaire","180°"]:GEO[X+"-"+(lp+rad1)];
  return{...base,lp,rad:rad1,exists:true,g,X,note:g?"":"Géométrie non gérée dans cette version."};}

function showMol(){const out=$("outMol"),m=analyse($("inMol").value);
  if(m.err){out.innerHTML=`<p class="bad">${m.err}</p>`;return;}
  let h=row("Molécule",subStr(m.formula))+row("Électrons de valence",`${m.total} au total · atome central : ${m.c} (${m.vc} e⁻ de valence)`);
  if(!m.exists){out.innerHTML=h+row("Existence",`<span class="bad">✗ N'existe pas</span><br>${m.why}`);return;}
  const ee=m.lp+(m.rad||0),notation=`AX${sub(m.X)}${ee?"E"+sub(ee):""}`;
  h+=row("Doublets",`${m.bo.reduce((a,b)=>a+b,0)} liants autour de ${m.c} · ${m.lp} non liant${m.lp>1?"s":""} sur ${m.c}`)
   +row("Notation VSEPR",notation);
  if(m.g)h+=row("Géométrie",`<b>${m.g[0]}</b> · angle ${m.g[1]}`);
  if(m.note)h+=row("Géométrie",m.note);
  h+=row("Schéma de Lewis","")+row("",drawMol(m))
   +row("Existence",m.rad?`<span class="ok">✓ Existe</span> mais c'est un radical (1 électron célibataire) : espèce très réactive.`:`<span class="ok">✓ Existe</span> : tous les atomes respectent leur capacité en électrons.`);
  out.innerHTML=h;}

/* ---------- Interface ---------- */
function tab(w){$("atome").hidden=w!=="a";$("molecule").hidden=w!=="m";
  $("tAtome").classList.toggle("on",w==="a");$("tMol").classList.toggle("on",w==="m");}
$("tAtome").onclick=()=>tab("a");$("tMol").onclick=()=>tab("m");
$("goAtome").onclick=showAtom;$("goMol").onclick=showMol;
$("inAtome").onkeydown=e=>e.key==="Enter"&&showAtom();
$("inCharge").onkeydown=e=>e.key==="Enter"&&showAtom();
$("inMol").onkeydown=e=>e.key==="Enter"&&showMol();
[["H",0],["C",0],["Fe",3],["Cu",0],["Cl",-1],["Na",1]].forEach(([s,q])=>{
  const b=document.createElement("button");b.textContent=s+(q?(q>0?"+":"−")+(Math.abs(q)>1?Math.abs(q):""):"");
  b.onclick=()=>{$("inAtome").value=s;$("inCharge").value=q;showAtom();};$("chipsA").append(b);});
"H2O,CO2,NH3,CH4,BF3,SO2,PCl5,SF6,XeF4,ClF3,NF5,CH5".split(",").forEach(f=>{
  const b=document.createElement("button");b.textContent=subStr(f);
  b.onclick=()=>{$("inMol").value=f;showMol();};$("chipsM").append(b);});