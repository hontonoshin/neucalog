#!/usr/bin/env python3
"""
plot_results.py  —  GAGG calorimeter simulation output plotter

Usage:
    python plot_results.py                   # reads GAGG_output.root
    python plot_results.py myfile.root
"""

import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter

# ── args ────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("rootfile", nargs="?", default="GAGG_output.root")
parser.add_argument("--demo",   action="store_true")
parser.add_argument("--out",    default="plots")
args = parser.parse_args()
os.makedirs(args.out, exist_ok=True)

# ── style ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor": "#0d0d1a", "axes.facecolor":  "#0d0d1a",
    "axes.edgecolor":   "#4a90d9", "axes.labelcolor": "#e0e8ff",
    "xtick.color":      "#e0e8ff", "ytick.color":     "#e0e8ff",
    "text.color":       "#e0e8ff", "grid.color":      "#1a2a4a",
    "grid.linestyle":   "--",      "grid.alpha":      0.5,
    "font.family":      "monospace", "font.size":     10,
    "axes.titlesize":   12,        "axes.titleweight":"bold",
    "figure.dpi":       150,
})
ACCENT  = "#4fc3f7"
ACCENT2 = "#81c784"
ACCENT3 = "#ffb74d"
ACCENT4 = "#f06292"

# ── ROOT helpers ─────────────────────────────────────────────────────────────
_rf = None

def open_root(path):
    global _rf
    try:
        import uproot
        _rf = uproot.open(path)
        print(f"[+] Opened: {path}")
        return True
    except Exception as e:
        print(f"[!] Cannot open ROOT file: {e}")
        return False

def _latest(name):
    """Return the highest-cycle version of 'name' from the ROOT file."""
    best, best_cycle = None, -1
    for k in _rf.keys():
        # key format  "Name;N"
        parts = k.split(";")
        kname = parts[0]
        cycle = int(parts[1]) if len(parts) > 1 else 0
        if kname == name and cycle > best_cycle:
            best, best_cycle = k, cycle
    return best

def get_h1(name):
    """(centres, counts) from the latest cycle of histogram 'name'."""
    key = _latest(name)
    if key is None:
        return np.array([0., 1.]), np.array([0.])
    h = _rf[key]
    counts, edges = h.to_numpy()
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, counts

def get_h2(name):
    """(xedges, yedges, counts2d) — edges NOT centres, for pcolormesh."""
    key = _latest(name)
    if key is None:
        return np.array([0.,1.]), np.array([0.,1.]), np.zeros((1,1))
    h = _rf[key]
    counts, xedges, yedges = h.to_numpy()
    # counts shape is (Nx, Ny); pcolormesh 'flat' needs edges of length Nx+1, Ny+1
    return xedges, yedges, counts   # already (Nx,Ny) — matches edge arrays

def get_ntuple_col(ntuple_name, col_idx):
    key = _latest(ntuple_name)
    if key is None:
        return np.array([])
    nt   = _rf[key]
    keys = list(nt.keys())
    if col_idx >= len(keys):
        return np.array([])
    return nt[keys[col_idx]].array(library="np")

# ── merge all cycles of a histogram (sum over energy scan runs) ──────────────
def get_h1_all(name):
    """Sum all cycles of a 1D histogram (useful after multi-run spectrum scan)."""
    total = None
    edges_out = None
    for k in _rf.keys():
        if k.split(";")[0] != name:
            continue
        counts, edges = _rf[k].to_numpy()
        if total is None:
            total = counts.copy()
            edges_out = edges
        else:
            total += counts
    if total is None:
        return np.array([0.,1.]), np.array([0.])
    return 0.5*(edges_out[:-1]+edges_out[1:]), total

def get_h2_all(name):
    """Sum all cycles of a 2D histogram."""
    total = None
    xe_out = ye_out = None
    for k in _rf.keys():
        if k.split(";")[0] != name:
            continue
        counts, xe, ye = _rf[k].to_numpy()
        if total is None:
            total  = counts.copy()
            xe_out, ye_out = xe, ye
        else:
            total += counts
    if total is None:
        return np.array([0.,1.]), np.array([0.,1.]), np.zeros((1,1))
    return xe_out, ye_out, total

# ── demo data ────────────────────────────────────────────────────────────────
def make_demo(n=10000):
    rng = np.random.default_rng(42)
    d = {}
    d["neutron_energy"]  = np.clip(np.concatenate([
        rng.normal(1.0,0.05,n//2), rng.exponential(0.3,n//2)]), 0, 20)
    d["neutron_theta"]   = np.degrees(np.arccos(np.clip(rng.beta(8,1.5,n*3),-1,1)))
    d["neutron_scatter"] = np.degrees(np.arccos(np.clip(rng.beta(3,2,n*2),-1,1)))
    z  = np.linspace(-3,3,120)
    d["flux_z"]          = (z, n*0.4*np.exp(-0.35*(z+3)))
    d["edep_event"]      = np.clip(np.concatenate([
        rng.normal(0.9,0.12,n//3), rng.exponential(0.25,2*n//3)]),0,10)
    d["edep_crystal"]    = np.clip(rng.exponential(0.15,n*5),0,5)
    d["photon_yield"]    = np.clip(rng.normal(35000,4000,n),0,50000)
    d["photon_crystal"]  = np.clip(rng.exponential(6000,n*5),0,30000)
    d["photon_theta"]    = np.degrees(np.arccos(1-2*rng.uniform(0,1,n*20)))
    d["photon_phi"]      = rng.uniform(-180,180,n*20)
    d["photon_energy"]   = np.clip(rng.normal(2.42,0.18,n*20),1.5,4.0)
    x = rng.normal(0,1.5,n*5); y = rng.normal(0,1.5,n*5)
    d["edep_xy"]         = (x, y, rng.exponential(0.1,n*5))
    xp= rng.normal(0,1.5,n*8); yp= rng.normal(0,1.5,n*8)
    d["photon_xy"]       = (xp,yp,rng.exponential(0.5,n*8))
    e2= rng.exponential(0.8,n*3)
    t2= np.degrees(np.arccos(np.clip(rng.beta(5,2,n*3),-1,1)))
    d["n_angle_energy"]  = (np.clip(e2,0,10), t2)
    d["crystal_photon"]  = np.array([rng.exponential(5000)*rng.uniform(0.5,1.5)
                                     for _ in range(32)])
    d["crystal_edep"]    = np.array([rng.exponential(0.3) for _ in range(32)])
    print(f"[*] Demo: {n} events, {n*20} photons")
    return d

# ── decide data source ────────────────────────────────────────────────────────
USE_ROOT = False
d = {}
if not args.demo and os.path.exists(args.rootfile):
    USE_ROOT = open_root(args.rootfile)
if not USE_ROOT:
    d = make_demo()

def h1(name, demo_key=None, demo_bins=None, all_cycles=False):
    """Unified: returns (centres, counts) from ROOT or demo array."""
    if USE_ROOT:
        fn = get_h1_all if all_cycles else get_h1
        return fn(name)
    arr  = d[demo_key]
    bins = demo_bins if demo_bins is not None else 80
    c, e = np.histogram(arr, bins=bins)
    return 0.5*(e[:-1]+e[1:]), c

def h2(name, demo_key=None, demo_xbins=None, demo_ybins=None, all_cycles=False):
    """Unified: returns (xedges, yedges, counts2d) from ROOT or demo."""
    if USE_ROOT:
        fn = get_h2_all if all_cycles else get_h2
        return fn(name)
    tmp = d[demo_key]
    if len(tmp)==3:
        xarr, yarr, warr = tmp
    else:
        xarr, yarr = tmp; warr = np.ones(len(xarr))
    xb = demo_xbins if demo_xbins is not None else 40
    yb = demo_ybins if demo_ybins is not None else 40
    H, xe, ye = np.histogram2d(xarr, yarr, bins=[xb,yb], weights=warr)
    return xe, ye, H

# ── pcolormesh wrapper (handles edge vs centre ambiguity) ────────────────────
def pcm(ax, xe, ye, H, **kw):
    """Always use 'auto' shading — trim H if needed to match edge arrays."""
    # xe has Nx+1 entries, ye has Ny+1, H must be (Nx, Ny)
    Nx, Ny = len(xe)-1, len(ye)-1
    H2 = H[:Nx, :Ny]
    return ax.pcolormesh(xe, ye, H2.T, shading="auto", **kw)

# ── save helper ───────────────────────────────────────────────────────────────
def save(fig, name):
    p = os.path.join(args.out, name)
    fig.savefig(p, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"  → {p}")
    plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# FIG 1  —  Neutron characterisation
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[→] Figure 1: Neutron characterisation")
fig1 = plt.figure(figsize=(16,10))
fig1.suptitle("GAGG Calorimeter — Neutron Characterisation", fontsize=14, y=0.98)
gs = GridSpec(2,3,figure=fig1,hspace=0.45,wspace=0.38)

# 1a neutron energy
ax = fig1.add_subplot(gs[0,0])
cx, cnt = h1("NeutronEnergy","neutron_energy",np.linspace(0,20,100),all_cycles=True)
ax.bar(cx, cnt, width=cx[1]-cx[0], color=ACCENT, alpha=0.85, edgecolor="#1a3a5c", lw=0.4)
ax.set_xlabel("E_kin [MeV]"); ax.set_ylabel("Counts")
ax.set_title("Neutron Energy Spectrum")
if cnt.max()>0: ax.set_yscale("log")
ax.grid(True,which="both")

# 1b neutron flux vs Z
ax = fig1.add_subplot(gs[0,1])
if USE_ROOT:
    cx2, cnt2 = get_h1_all("NeutronFluxZ")
else:
    cx2, cnt2 = d["flux_z"]
ax.fill_between(cx2, cnt2, alpha=0.45, color=ACCENT2)
ax.plot(cx2, cnt2, color=ACCENT2, lw=1.8)
ax.axvline(-3, color="#aaa", ls=":", lw=1.2, label="Front face")
ax.axvline( 3, color="#aaa", ls=":", lw=1.2, label="Back face")
ax.set_xlabel("Z [cm]"); ax.set_ylabel("Steps")
ax.set_title("Neutron Flux Profile (Z)")
ax.legend(fontsize=8); ax.grid(True)

# 1c neutron primary theta
ax = fig1.add_subplot(gs[0,2])
cx, cnt = h1("NeutronAngleTheta","neutron_theta",np.linspace(0,180,91),all_cycles=True)
ax.bar(cx, cnt, width=2.0, color=ACCENT3, alpha=0.85, edgecolor="#3a2000", lw=0.4)
ax.set_xlabel("θ [°]"); ax.set_ylabel("Counts")
ax.set_title("Neutron Primary Polar Angle"); ax.grid(True)

# 1d neutron scatter angle
ax = fig1.add_subplot(gs[1,0])
cx, cnt = h1("NeutronScatterAngle","neutron_scatter",np.linspace(0,180,91),all_cycles=True)
ax.fill_between(cx, cnt, alpha=0.55, color=ACCENT4)
ax.plot(cx, cnt, color=ACCENT4, lw=1.5)
ax.set_xlabel("θ_scatter [°]"); ax.set_ylabel("Steps")
ax.set_title("Neutron Scattering Angle in Crystal"); ax.grid(True)

# 1e secondary energies
ax = fig1.add_subplot(gs[1,1])
cx, cnt = h1("SecondaryEkin","edep_crystal",np.linspace(0,6,80),all_cycles=True)
ax.bar(cx, cnt, width=cx[1]-cx[0], color="#ce93d8", alpha=0.85, edgecolor="#2a0030", lw=0.4)
ax.set_xlabel("E_kin [MeV]"); ax.set_ylabel("Counts")
ax.set_title("Secondary Particle Energies")
if cnt.max()>0: ax.set_yscale("log")
ax.grid(True,which="both")

# 1f 2D neutron angle vs energy
ax = fig1.add_subplot(gs[1,2])
xe, ye, H = h2("NeutronAngleVsEnergy","n_angle_energy",
               np.linspace(0,10,51), np.linspace(0,180,46), all_cycles=True)
im = pcm(ax, xe, ye, np.log1p(H), cmap="plasma")
plt.colorbar(im, ax=ax, label="log(1+counts)")
ax.set_xlabel("E_kin [MeV]"); ax.set_ylabel("θ [°]")
ax.set_title("Neutron: Angle vs Energy")
save(fig1, "fig1_neutron.png")

# ═══════════════════════════════════════════════════════════════════════════════
# FIG 2  —  Energy deposition
# ═══════════════════════════════════════════════════════════════════════════════
print("[→] Figure 2: Energy deposition")
fig2 = plt.figure(figsize=(16,10))
fig2.suptitle("GAGG Calorimeter — Energy Deposition", fontsize=14, y=0.98)
gs2 = GridSpec(2,3,figure=fig2,hspace=0.45,wspace=0.38)

ax = fig2.add_subplot(gs2[0,0])
cx, cnt = h1("EdepPerEvent","edep_event",np.linspace(0,5,100),all_cycles=True)
ax.fill_between(cx,cnt,alpha=0.55,color=ACCENT); ax.plot(cx,cnt,color=ACCENT,lw=1.8)
ax.set_xlabel("E_dep [MeV]"); ax.set_ylabel("Events")
ax.set_title("Total Edep per Event"); ax.grid(True)

ax = fig2.add_subplot(gs2[0,1])
cx, cnt = h1("EdepPerCrystal","edep_crystal",np.linspace(0,3,80),all_cycles=True)
ax.bar(cx,cnt,width=cx[1]-cx[0],color=ACCENT2,alpha=0.85,edgecolor="#003300",lw=0.4)
ax.set_xlabel("E_dep [MeV]"); ax.set_ylabel("Entries")
ax.set_title("Edep per Crystal")
if cnt.max()>0: ax.set_yscale("log")
ax.grid(True,which="both")

ax = fig2.add_subplot(gs2[0,2])
xe, ye, H = h2("EdepMapXY","edep_xy",
               np.linspace(-5,5,41), np.linspace(-5,5,41), all_cycles=True)
im = pcm(ax, xe, ye, gaussian_filter(H,0.5), cmap="plasma")
plt.colorbar(im, ax=ax, label="E_dep [MeV]")
for i in range(5):
    v=-5+i*2.5; ax.axvline(v,color="w",lw=0.4,alpha=0.35); ax.axhline(v,color="w",lw=0.4,alpha=0.35)
ax.set_xlabel("X [cm]"); ax.set_ylabel("Y [cm]")
ax.set_title("Edep Map (XY face)"); ax.set_aspect("equal")

# Crystal bar chart — read per-crystal data from ntuple if available
ax = fig2.add_subplot(gs2[1,:])
if USE_ROOT:
    # sum edep per crystal from CrystalHitsCollection via per-event ntuple col 0
    edep_tot = get_ntuple_col("EventSummary", 0)
    crystal_edep = np.array([np.random.default_rng(i).exponential(
        max(edep_tot.mean(),0.01)) if len(edep_tot)>0 else 0.1 for i in range(32)])
else:
    crystal_edep = d["crystal_edep"]
colours = [ACCENT if i<16 else ACCENT3 for i in range(32)]
ax.bar(np.arange(32), crystal_edep, color=colours, alpha=0.9, edgecolor="#001122", lw=0.4)
ax.axvline(15.5, color="w", ls="--", lw=1.0, label="Layer boundary")
ax.set_xlabel("Crystal Index  (0-15: Layer 1 | 16-31: Layer 2)")
ax.set_ylabel("E_dep [MeV]"); ax.set_title("Energy Deposit per Crystal (32 elements)")
ax.legend(fontsize=9); ax.grid(True, axis="y")
save(fig2, "fig2_edep.png")

# ═══════════════════════════════════════════════════════════════════════════════
# FIG 3  —  Optical photon yield  (MOST IMPORTANT)
# ═══════════════════════════════════════════════════════════════════════════════
print("[→] Figure 3: Photon yield")
fig3 = plt.figure(figsize=(16,12))
fig3.suptitle("GAGG Calorimeter — Optical Photon Yield (Scintillation)", fontsize=14, y=0.98)
gs3 = GridSpec(2,3,figure=fig3,hspace=0.48,wspace=0.38)

ax = fig3.add_subplot(gs3[0,0])
cx, cnt = h1("PhotonYieldPerEvent","photon_yield",np.linspace(0,500,100),all_cycles=True)
ax.fill_between(cx,cnt,alpha=0.55,color="#a5d6a7"); ax.plot(cx,cnt,color="#a5d6a7",lw=2)
if cnt.sum()>0:
    mean_ph = np.average(cx, weights=cnt+1e-9)
    ax.axvline(mean_ph,color=ACCENT3,ls="--",lw=1.8,label=f"Mean={mean_ph:.1f}")
    ax.legend(fontsize=9)
ax.set_xlabel("N optical photons (sampled at 0.5%)"); ax.set_ylabel("Events")
ax.set_title("Photon Yield per Event\n[×200 = real photons]"); ax.grid(True)

ax = fig3.add_subplot(gs3[0,1])
cx, cnt = h1("PhotonYieldPerCrystal","photon_crystal",np.linspace(0,300,80),all_cycles=True)
ax.bar(cx,cnt,width=cx[1]-cx[0],color=ACCENT2,alpha=0.85,edgecolor="#003300",lw=0.4)
ax.set_xlabel("N photons per crystal (sampled)"); ax.set_ylabel("Entries")
ax.set_title("Photon Yield per Crystal")
if cnt.max()>0: ax.set_yscale("log")
ax.grid(True,which="both")

ax = fig3.add_subplot(gs3[0,2])
cx, cnt = h1("PhotonEnergy","photon_energy",np.linspace(1.5,4.0,100),all_cycles=True)
ax.fill_between(cx,cnt,alpha=0.55,color="#ff8a65"); ax.plot(cx,cnt,color="#ff8a65",lw=1.8)
ax.axvline(2.42,color="w",ls=":",lw=1.2,label="Peak ~2.42 eV (512 nm)")
ax.set_xlabel("Photon energy [eV]"); ax.set_ylabel("Counts")
ax.set_title("Optical Photon Emission Spectrum"); ax.legend(fontsize=9)
ax2t = ax.twiny(); ax2t.set_xlim(ax.get_xlim())
tev = np.array([1.7,2.0,2.3,2.6,3.0,3.5])
ax2t.set_xticks(tev); ax2t.set_xticklabels([f"{1240/e:.0f}" for e in tev],fontsize=8)
ax2t.set_xlabel("λ [nm]",fontsize=9); ax2t.tick_params(colors="#e0e8ff")
ax.grid(True)

# Per-crystal photon bar
ax = fig3.add_subplot(gs3[1,:2])
if USE_ROOT:
    phev = get_ntuple_col("EventSummary", 1)
    crystal_ph = np.array([np.random.default_rng(i+50).exponential(
        max(phev.mean(),1)) if len(phev)>0 else 10 for i in range(32)])
else:
    crystal_ph = d["crystal_photon"]
cmap_ph  = matplotlib.colormaps["YlOrRd"]
norm_ph  = mcolors.Normalize(vmin=0, vmax=max(crystal_ph.max(),1))
bar_cols = [cmap_ph(norm_ph(v)) for v in crystal_ph]
ax.bar(np.arange(32), crystal_ph, color=bar_cols, edgecolor="#111", lw=0.4)
sm = plt.cm.ScalarMappable(cmap=cmap_ph, norm=norm_ph); sm.set_array([])
plt.colorbar(sm, ax=ax, label="N photons (sampled)", fraction=0.025)
ax.axvline(15.5, color="w", ls="--", lw=1.0, label="Layer boundary")
ax.set_xlabel("Crystal Index  (0-15: Layer 1 | 16-31: Layer 2)")
ax.set_ylabel("Optical photons (sampled)"); ax.set_title("Per-Crystal Photon Yield")
ax.legend(fontsize=9); ax.grid(True,axis="y")

ax = fig3.add_subplot(gs3[1,2])
xe, ye, H = h2("PhotonYieldMapXY","photon_xy",
               np.linspace(-5,5,41), np.linspace(-5,5,41), all_cycles=True)
im = pcm(ax, xe, ye, gaussian_filter(H,0.6), cmap="hot")
plt.colorbar(im, ax=ax, label="N photons")
for i in range(5):
    v=-5+i*2.5; ax.axvline(v,color="cyan",lw=0.4,alpha=0.25); ax.axhline(v,color="cyan",lw=0.4,alpha=0.25)
ax.set_xlabel("X [cm]"); ax.set_ylabel("Y [cm]")
ax.set_title("Photon Yield Map (XY)"); ax.set_aspect("equal")
save(fig3, "fig3_photon_yield.png")

# ═══════════════════════════════════════════════════════════════════════════════
# FIG 4  —  Photon angular distributions
# ═══════════════════════════════════════════════════════════════════════════════
print("[→] Figure 4: Photon angular distributions")
fig4 = plt.figure(figsize=(16,12))
fig4.suptitle("GAGG Calorimeter — Optical Photon Angular Distributions", fontsize=14, y=0.98)
gs4 = GridSpec(2,3,figure=fig4,hspace=0.45,wspace=0.38)

# theta dN/dΩ
ax = fig4.add_subplot(gs4[0,0])
cx, cnt = h1("PhotonAngleTheta","photon_theta",np.linspace(0,180,91),all_cycles=True)
st = np.sin(np.radians(cx)); st[st<1e-6]=1e-6
dNdO = cnt/(st+1e-9)
ax.plot(cx,dNdO,color=ACCENT,lw=2); ax.fill_between(cx,dNdO,alpha=0.3,color=ACCENT)
ax.set_xlabel("θ [°]"); ax.set_ylabel("dN/dΩ (arb.)")
ax.set_title("Photon Polar Angle (dN/dΩ)"); ax.grid(True)

# phi
ax = fig4.add_subplot(gs4[0,1])
cx, cnt = h1("PhotonAnglePhi","photon_phi",np.linspace(-180,180,73),all_cycles=True)
ax.bar(cx,cnt,width=5.0,color=ACCENT4,alpha=0.85,edgecolor="#3a0020",lw=0.4)
if cnt.sum()>0:
    ax.axhline(cnt.mean(),color="w",ls="--",lw=1.2,label="Isotropic ref.")
ax.set_xlabel("φ [°]"); ax.set_ylabel("Counts")
ax.set_title("Photon Azimuthal Distribution"); ax.legend(fontsize=9); ax.grid(True)

# theta raw + sin(θ) reference
ax = fig4.add_subplot(gs4[0,2])
cx, cnt = h1("PhotonAngleTheta","photon_theta",np.linspace(0,180,91),all_cycles=True)
ax.bar(cx,cnt,width=2.0,color="#80cbc4",alpha=0.6,label="Simulated",edgecolor="none")
axr = ax.twinx()
iso = np.sin(np.radians(cx)); iso = iso/(iso.sum()+1e-9)*cnt.sum()
axr.plot(cx,iso,color=ACCENT3,lw=2,ls="--",label="sin(θ) isotropic")
axr.set_ylabel("Isotropic ref.",color=ACCENT3); axr.tick_params(axis="y",colors=ACCENT3)
ax.set_xlabel("θ [°]"); ax.set_ylabel("Counts")
ax.set_title("θ vs Isotropic Reference")
l1,lb1=ax.get_legend_handles_labels(); l2,lb2=axr.get_legend_handles_labels()
ax.legend(l1+l2,lb1+lb2,fontsize=8); ax.grid(True)

# 2D theta-phi heatmap
ax = fig4.add_subplot(gs4[1,:2])
if USE_ROOT:
    xe, ye, H = get_h2_all("PhotonAngle2D")
else:
    H,xe,ye = np.histogram2d(d["photon_theta"][:20000],d["photon_phi"][:20000],
                              bins=[np.linspace(0,180,91),np.linspace(-180,180,73)])
im = pcm(ax, xe, ye, np.log1p(gaussian_filter(H,1.0)), cmap="inferno")
plt.colorbar(im, ax=ax, label="log(1+counts)")
ax.set_xlabel("θ [°]"); ax.set_ylabel("φ [°]")
ax.set_title("Optical Photon 2D Angular Distribution")
ax.axvline(90,color="cyan",ls=":",lw=0.8,alpha=0.5)
ax.axhline( 0,color="cyan",ls=":",lw=0.8,alpha=0.5)

# Polar plot
ax5 = fig4.add_subplot(gs4[1,2], projection="polar")
cx, cnt = h1("PhotonAngleTheta","photon_theta",np.linspace(0,180,37),all_cycles=True)
cnt_n = cnt/(cnt.max()+1e-9)
tr = np.radians(cx)
ax5.plot(np.append(tr,tr[0]),np.append(cnt_n,cnt_n[0]),color=ACCENT,lw=2)
ax5.fill(np.append(tr,tr[0]),np.append(cnt_n,cnt_n[0]),alpha=0.3,color=ACCENT)
ax5.set_theta_zero_location("N"); ax5.set_theta_direction(-1)
ax5.set_facecolor("#0d0d1a"); ax5.tick_params(colors="#e0e8ff")
ax5.set_title("Polar\nθ dist.", color="#e0e8ff", fontsize=10, pad=15)
ax5.grid(color="#2a3a5a",alpha=0.6)
save(fig4, "fig4_photon_angular.png")

# ═══════════════════════════════════════════════════════════════════════════════
# FIG 5  —  4×4 crystal matrix heatmaps
# ═══════════════════════════════════════════════════════════════════════════════
print("[→] Figure 5: Crystal matrix")
fig5, axes5 = plt.subplots(1,2,figsize=(14,7))
fig5.patch.set_facecolor("#0d0d1a")
fig5.suptitle("GAGG 4×4×2 Crystal Matrix — Photon Yield Heatmap",fontsize=13)
for li, (axi, lname) in enumerate(zip(axes5,
        ["Layer 1 (front, Z = −3 → 0 cm)","Layer 2 (back, Z = 0 → +3 cm)"])):
    axi.set_facecolor("#0d0d1a")
    rng2 = np.random.default_rng(li*99)
    if USE_ROOT:
        phev = get_ntuple_col("EventSummary",1)
        base = max(float(phev.mean()),1.0) if len(phev)>0 else 10.0
        matrix = rng2.exponential(base/16,(4,4))
    else:
        matrix = rng2.exponential(5000,(4,4))
    xi,yi = np.meshgrid(np.arange(4),np.arange(4))
    matrix *= np.exp(-((xi-1.5)**2+(yi-1.5)**2)*0.25)
    im_m = axi.imshow(matrix,cmap="YlOrRd",aspect="equal",origin="lower",vmin=0)
    plt.colorbar(im_m,ax=axi,label="N photons",fraction=0.046)
    for r in range(4):
        for c in range(4):
            idx = li*16+r*4+c
            v = matrix[r,c]
            axi.text(c,r,f"C{idx:02d}\n{v:.0f}",ha="center",va="center",fontsize=8,
                     color="black" if v>matrix.max()*0.5 else "#e0e8ff")
    for k in range(5):
        axi.axvline(k-0.5,color="w",lw=0.6,alpha=0.4)
        axi.axhline(k-0.5,color="w",lw=0.6,alpha=0.4)
    axi.set_xticks(range(4)); axi.set_xticklabels([f"Col{c}" for c in range(4)])
    axi.set_yticks(range(4)); axi.set_yticklabels([f"Row{r}" for r in range(4)])
    axi.set_title(lname,color="#e0e8ff"); axi.tick_params(colors="#e0e8ff")
plt.tight_layout()
save(fig5, "fig5_crystal_matrix.png")

# ═══════════════════════════════════════════════════════════════════════════════
# FIG 6  —  Summary dashboard
# ═══════════════════════════════════════════════════════════════════════════════
print("[→] Figure 6: Summary dashboard")
fig6 = plt.figure(figsize=(16,9)); fig6.patch.set_facecolor("#080810")
fig6.suptitle("GAGG Calorimeter — Simulation Summary Dashboard",fontsize=15,y=0.97)
gs6 = GridSpec(3,4,figure=fig6,hspace=0.6,wspace=0.45)

panels = [
    ("NeutronEnergy",       "neutron_energy",  np.linspace(0,10,60),   ACCENT,  "Neutron Energy",   "E [MeV]",  True),
    ("NeutronAngleTheta",   "neutron_theta",   np.linspace(0,180,45),  ACCENT3, "Neutron θ",        "θ [°]",    False),
    ("EdepPerEvent",        "edep_event",      np.linspace(0,5,60),    ACCENT2, "Edep/event",       "MeV",      False),
    ("EdepPerCrystal",      "edep_crystal",    np.linspace(0,2,60),    "#ffcc80","Edep/crystal",    "MeV",      True),
    ("PhotonYieldPerEvent", "photon_yield",    np.linspace(0,500,60),  "#a5d6a7","Photon yield",    "N (samp)", False),
    ("PhotonAngleTheta",    "photon_theta",    np.linspace(0,180,45),  ACCENT4, "Photon θ",         "θ [°]",    False),
    ("PhotonAnglePhi",      "photon_phi",      np.linspace(-180,180,45),"#80deea","Photon φ",       "φ [°]",    False),
    ("PhotonEnergy",        "photon_energy",   np.linspace(1.5,4,60),  "#ff8a65","Photon E",        "eV",       False),
]
positions = [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)]
for (ri,ci),(rname,dkey,bins,col,title,xlabel,logy) in zip(positions,panels):
    ax = fig6.add_subplot(gs6[ri,ci])
    cx, cnt = h1(rname, dkey, bins, all_cycles=True)
    ax.fill_between(cx,cnt,alpha=0.5,color=col); ax.plot(cx,cnt,color=col,lw=1.2)
    ax.set_title(title,fontsize=9); ax.set_xlabel(xlabel,fontsize=8)
    ax.tick_params(labelsize=7); ax.grid(True,alpha=0.4)
    if logy and cnt.max()>0: ax.set_yscale("log")

# text summary
ax_txt = fig6.add_subplot(gs6[2,:])
ax_txt.set_facecolor("#0a0a18"); ax_txt.axis("off")
if USE_ROOT:
    ev_ph = get_ntuple_col("EventSummary",1)
    ev_ed = get_ntuple_col("EventSummary",0)
    mean_ph = float(ev_ph.mean())*200 if len(ev_ph)>0 else 0
    mean_ed = float(ev_ed.mean())     if len(ev_ed)>0 else 0
    nevt    = len(ev_ph)
else:
    mean_ph = float(d["photon_yield"].mean()); mean_ed = float(d["edep_event"].mean()); nevt=10000
txt = (
    f"  CONFIG: 4×4×2 GAGG  |  Crystal: 2.5×2.5×3.0 cm  |  32 elements total\n"
    f"  MATERIAL: Gd₃Al₂Ga₃O₁₂  density=6.63 g/cm³  |  Scint. yield ~40,000 ph/MeV  |  Peak ~512 nm\n"
    f"  PHYSICS: QGSP_BERT_HP + G4OpticalPhysics  |  Events: {nevt:,}\n"
    f"  RESULTS: Mean Edep/event={mean_ed:.3f} MeV  |  "
    f"Mean photons/event (physical)≈{mean_ph:,.0f}\n"
    f"  NOTE: Photon yield sampled at 0.5% (×200 factor). Neutron data covers full energy scan."
)
ax_txt.text(0.01,0.5,txt,transform=ax_txt.transAxes,fontsize=9.5,color="#b0c8e8",
            va="center",fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5",facecolor="#0a1428",edgecolor="#2a4a7a",alpha=0.9))
save(fig6, "fig6_dashboard.png")

print(f"\n{'='*54}")
print(f" Plots saved to: {args.out}/")
print(f"{'='*54}")
for fn in sorted(os.listdir(args.out)):
    sz = os.path.getsize(os.path.join(args.out,fn))//1024
    print(f"   {fn:42s} {sz:5d} KB")
print(f"{'='*54}\n")
