"""Calculos propios del capitulo 4 (termica en vacio y blindaje)."""
import numpy as np
c=299792458.0; mu0=4e-7*np.pi; sb=5.670374419e-8
kcu=387.6; kfr=0.35; kal=218.0           # Guerra et al.
print("== PCB conductividad efectiva (1,6 mm)")
for ncu,tcu_um in [(4,35),(4,70),(2,35)]:
    t=1.6e-3; tc=ncu*tcu_um*1e-6; tf=t-tc
    kpar=(kcu*tc+kfr*tf)/t; kperp=t/(tc/kcu+tf/kfr)
    print(f"{ncu} capas Cu {tcu_um}um: Cu total={tc*1e3:.3f}mm  k_par={kpar:.1f}  k_perp={kperp:.3f}")
print("== vias")
def Rvia(d_mm,tp_um,L_mm=1.6):
    ro=d_mm/2*1e-3; ri=ro-tp_um*1e-6; A=np.pi*(ro**2-ri**2); return L_mm*1e-3/(kcu*A), A
for tp in (25,35):
    R,A=Rvia(0.3,tp); print(f"via 0.3mm tp={tp}um A={A*1e6:.4f}mm2 R={R:.0f} K/W; x9={R/9:.1f}; x25={R/25:.1f}")
Rfr=1.6e-3/(kfr*25e-6); print("FR-4 5x5mm slab R=",round(Rfr,1))
R25,_=Rvia(0.3,25)
for n in (1,5,9,16,25): print(n, round(1/(n/R25+1/Rfr),1))
print("== radiacion vs conduccion")
A=0.0902*0.0959; 
for eps in (0.9,0.5,0.1):
    hr=4*eps*sb*300**3; print(f"eps={eps} h_rad={hr:.2f} W/m2K  G(2 caras)={2*A*hr:.3f} W/K")
od,id_,L=5e-3,3.2e-3,15e-3; As=np.pi/4*(od**2-id_**2); Gs=kal*As/L; print("spacer A=",As*1e6,"mm2 G=",round(Gs,3),"W/K")
for Gj in (0.105,0.21):
    R=1/Gs+2/Gj; G1=1/R; print(f"junta {Gj}: R_path={R:.1f} K/W  G4={4*G1:.3f} W/K  dT@1W={1/(4*G1):.2f} K  dT@2W={2/(4*G1):.2f}")
print("== blindaje")
def delta(f,s): return 1/np.sqrt(np.pi*f*mu0*s)
for nm,s in (("Cu",5.8e7),("Al",3.7e7)):
    for f in (145e6,435e6):
        d=delta(f,s); Zs=np.sqrt(2*np.pi*f*mu0/s)
        for t in (0.2e-3,):
            Aab=8.686*t/d; R=20*np.log10(376.73/(4*Zs))
            print(f"{nm} f={f/1e6:.0f} delta={d*1e6:.2f}um Zs={Zs*1e3:.2f}mOhm A(0.2mm)={Aab:.0f}dB R_plane={R:.1f}dB")
print("== aperturas")
for l in (1,2,5,10,20):
    fc=15000/(l/10)  # MHz (l en cm)
    print(f"l={l}mm fc={fc/1000:.1f}GHz  SE@145={20*np.log10(fc/145):.1f}  SE@435={20*np.log10(fc/435):.1f}  SE@870={20*np.log10(fc/870):.1f}")
print("lambda/20:", [round(c/f/20*1e3,1) for f in (435e6,870e6,1305e6)])
print("SE a paso lambda/20 =",20*np.log10(10))
for p in (34,17,10,5,2): print(p,"mm",[round(20*np.log10(c/f*1e3/(2*p)),1) for f in (435e6,870e6)])
print("N huecos 10log:", [round(10*np.log10(n),1) for n in (2,4,10)])
print("ab hueco 27.3t/l t=0.2mm l=5mm:",27.3*0.2/5)

# ---------------- Figuras ----------------
import os, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
out=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','Figuras')
com=FuncFormatter(lambda v,p: f"{v:g}".replace('.',','))
l=np.linspace(1,40,200)
fig,ax=plt.subplots(figsize=(6,3.3))
for f in (145,435,870):
    ax.plot(l,np.clip(20*np.log10(15000/(l/10)/f),0,None),label=f"{f} MHz")
ax.set_xlabel("Dimensión mayor de la abertura [mm]"); ax.set_ylabel("Efectividad de blindaje [dB]")
ax.xaxis.set_major_formatter(com); ax.yaxis.set_major_formatter(com); ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(out,'se_apertura.pdf')); plt.close(fig)
n=np.arange(1,31)
fig,ax=plt.subplots(figsize=(6,3.3))
for tp,col in ((25,'C0'),(35,'C1')):
    R,_=Rvia(0.3,tp); ax.plot(n,1/(n/R+1/Rfr),color=col,label=f"vías de 0,3 mm, metalizado de {tp} µm")
ax.axhline(Rfr,color='k',ls='--',lw=1,label="solo FR-4 bajo el pad de 5 mm × 5 mm")
ax.set_xlabel("Cantidad de vías bajo el pad"); ax.set_ylabel("Resistencia térmica [K/W]")
ax.xaxis.set_major_formatter(com); ax.yaxis.set_major_formatter(com); ax.grid(alpha=.3); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(out,'resistencia_vias.pdf')); plt.close(fig)
print("figuras ok")
