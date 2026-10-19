"""Calculos propios del capitulo 3 (layout de PCB para RF).

Reproduce las tablas y figuras del capitulo:
  - microstrip de 50 ohm (Hammerstad-Jensen segun Steer, LibreTexts 4.4)
  - perdidas dielectricas y de conductor (formulas habituales, ver Pozar sec. 3.8)
  - CPWG (expresiones de Ghione y Naldi, tal como las implementa KiCad)
  - inductancia de via (aproximacion de Johnson)
Uso:  python calculos_rf.py        (requiere numpy, scipy y matplotlib)
Genera ../Figuras/perdidas_microstrip.pdf y ../Figuras/reactancia_vias.pdf
"""
import numpy as np
from scipy.special import ellipk
from scipy.optimize import brentq
c=299792458.0; mu0=4e-7*np.pi; eta=376.730313
def hj(u,er):  # Steer / Hammerstad-Jensen
    a=1+np.log((u**4+(u/52)**2)/(u**4+0.432))/49+np.log(1+(u/18.1)**3)/18.7
    b=0.564*((er-0.9)/(er+3))**0.053
    ee=(er+1)/2+(er-1)/2*(1+10/u)**(-a*b)
    F1=6+(2*np.pi-6)*np.exp(-(30.666/u)**0.7528)
    Z01=60*np.log(F1/u+np.sqrt(1+(2/u)**2))
    return ee,Z01/np.sqrt(ee)
# check Steer example: w=600,h=635,er=4.1 -> ee 2.967, Z0 75.4
print("check",hj(600/635,4.1))
def wh50(er,Z=50.0): return brentq(lambda u:hj(u,er)[1]-Z,0.05,20)
print("== microstrip 50 ohm")
for name,er,hs in [("FR4 1.6",4.0,1.6),("FR4 0.8",4.0,0.8),("FR4 0.2",4.0,0.2),("RO4350B 0.508",3.66,0.508),("RO4350B 0.762",3.66,0.762)]:
    u=wh50(er); ee,Z=hj(u,er); print(f"{name}: w/h={u:.3f} w={u*hs:.3f} mm ee={ee:.3f} Z={Z:.2f}")
# lengths
print("== lengths")
for er,hs in [(4.0,1.6),(3.66,0.508)]:
    u=wh50(er); ee,_=hj(u,er)
    for f in (145e6,435e6):
        l0=c/f*1e3; lg=l0/np.sqrt(ee)
        print(f"er={er} h={hs} f={f/1e6:.0f}MHz l0={l0:.0f}mm lg={lg:.0f}mm lg/20={lg/20:.1f} lg/10={lg/10:.1f} 90.2mm={90.2/lg:.3f}lg 95.9={95.9/lg:.3f}lg")
# losses
def rs(f,sig=5.8e7): return np.sqrt(np.pi*f*mu0/sig)
def delta(f,sig=5.8e7): return 1/np.sqrt(np.pi*f*mu0*sig)
print("== loss")
for name,er,hs,td in [("FR4 1.6",4.0,1.6,0.0161),("RO4350B 0.508",3.66,0.508,0.0037)]:
    u=wh50(er); ee,Z=hj(u,er); w=u*hs*1e-3
    for f in (145e6,435e6):
        k0=2*np.pi*f/c
        ad=k0*er*(ee-1)*td/(2*np.sqrt(ee)*(er-1))      # Np/m
        ac=rs(f)/(Z*w)
        print(f"{name} {f/1e6:.0f}MHz delta={delta(f)*1e6:.2f}um Rs={rs(f)*1e3:.2f}mohm ad={ad*8.686:.3f}dB/m ac={ac*8.686:.3f}dB/m total/10cm={(ad+ac)*8.686*0.1:.4f}dB")
# CPWG
def cpwg(W,S,H,er):
    k1=W/(W+2*S); q1=ellipk(k1**2)/ellipk(1-k1**2)
    k3=np.tanh(np.pi/4*W/H)/np.tanh(np.pi/4*(W+2*S)/H); q3=ellipk(k3**2)/ellipk(1-k3**2)
    qz=1/(q1+q3); ee=1+q3*qz*(er-1); Z=(eta/2)*qz/np.sqrt(ee); return ee,Z
print("== CPWG (zero thickness)")
for er,H,W in [(4.0,1.6,1.0),(4.0,1.6,0.5),(4.0,0.8,0.5),(3.66,0.508,0.5)]:
    S=brentq(lambda S:cpwg(W,S,H,er)[1]-50,0.02,5)
    ee,Z=cpwg(W,S,H,er); print(f"er={er} H={H} W={W} -> S={S:.3f} mm ee={ee:.3f} Z={Z:.2f}")
# CPWG limit check: large S -> microstrip
print("limit", cpwg(3.0,50,1.6,4.0), hj(3.0/1.6,4.0))
# vias
print("== vias")
def Lvia(h_mm,d_mm):
    h=h_mm/25.4; d=d_mm/25.4; return 5.08*h*(np.log(4*h/d)+1)
for h,d in [(1.6,0.3),(0.8,0.3),(1.6,0.2)]:
    L=Lvia(h,d)
    for f in (145e6,435e6): print(f"h={h} d={d} L={L:.2f}nH X@{f/1e6:.0f}={2*np.pi*f*L*1e-9:.2f}ohm, 4 vias X={2*np.pi*f*L/4*1e-9:.2f}")
print("== CPWG thin prepreg")
for er,H,W in [(4.0,0.2,0.4),(4.0,0.2,0.3),(4.0,0.2,0.25),(4.0,0.1,0.2)]:
    S=brentq(lambda S:cpwg(W,S,H,er)[1]-50,0.01,5)
    ee,Z=cpwg(W,S,H,er); print(f"er={er} H={H} W={W} -> S={S:.3f} mm ee={ee:.3f}")
print("ms 0.2mm w=",wh50(4.0)*0.2)
print("cpw sens 1.6/1.0 gap +-0.02:", [round(cpwg(1.0,s,1.6,4.0)[1],1) for s in (0.10,0.12,0.14)])
print("ms 1.6 w+-0.05:", [round(hj((3.289+d)/1.6,4.0)[1],2) for d in (-0.05,0,0.05)])
print("er sens ms 1.6:", [round(hj(3.289/1.6,e)[1],2) for e in (3.9,4.0,4.4,4.6)])
print("== CPWG 1.6 wide")
for W in (1.5,2.0,2.5):
    S=brentq(lambda S:cpwg(W,S,1.6,4.0)[1]-50,0.02,8); print(W,round(S,3),round(cpwg(W,S,1.6,4.0)[0],3),[round(cpwg(W,S+d,1.6,4.0)[1],1) for d in (-0.05,0.05)])


# ---------------- Figuras ----------------
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os
out=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','Figuras')

com=FuncFormatter(lambda v,p: f"{v:g}".replace('.',','))
def style(ax): 
    ax.xaxis.set_major_formatter(com); ax.yaxis.set_major_formatter(com); ax.grid(alpha=.3)
f=np.linspace(100e6,1000e6,200)
def loss(er,hs,td):
    u=brentq(lambda u:hj(u,er)[1]-50,0.05,20); ee,Z=hj(u,er); w=u*hs*1e-3
    k0=2*np.pi*f/c; ad=k0*er*(ee-1)*td/(2*np.sqrt(ee)*(er-1)); ac=np.sqrt(np.pi*f*mu0/5.8e7)/(Z*w)
    return ad*8.686*0.1, ac*8.686*0.1
fig,ax=plt.subplots(figsize=(6,3.4))
for lab,er,hs,td,col in [("FR-4 (FR406), 1,6 mm",4.0,1.6,0.0161,'C0'),("RO4350B, 0,508 mm",3.66,0.508,0.0037,'C3')]:
    ad,ac=loss(er,hs,td); ax.plot(f/1e6,ad+ac,color=col,label=lab+" (total)"); ax.plot(f/1e6,ad,color=col,ls='--',lw=1,label=lab+" (dieléctrica)")
ax.set_xlabel("Frecuencia [MHz]"); ax.set_ylabel("Atenuación en 10 cm [dB]"); style(ax); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(out+'/perdidas_microstrip.pdf'); plt.close(fig)
fig,ax=plt.subplots(figsize=(6,3.2))
L=5.08*(1.6/25.4)*(np.log(4*(1.6/25.4)/(0.3/25.4))+1)
for n in (1,2,4): ax.plot(f/1e6,2*np.pi*f*L*1e-9/n,label=f"{n} vía"+("s" if n>1 else "")+" en paralelo")
ax.set_xlabel("Frecuencia [MHz]"); ax.set_ylabel("Reactancia [Ω]"); style(ax); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(out+'/reactancia_vias.pdf'); plt.close(fig)
print("ok",L)
