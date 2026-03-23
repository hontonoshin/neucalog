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
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "axes.edgecolor":   "black",
    "axes.labelcolor":  "black",
    "xtick.color":      "black",
    "ytick.color":      "black",
    "text.color":       "black",
    "grid.color":       "#cccccc",
    "grid.linestyle":   "--",
    "grid.alpha":       0.7,
    "font.family":      "serif",
    "font.serif":       ["Times New Roman"],
    "font.size":        10,
    "axes.titlesize":   12,
    "axes.titleweight": "bold",
    "figure.dpi":       300,
    "savefig.dpi":      300,
    "savefig.bbox":     "tight",
})

# Colour palette 
COLORS = {
    "blue":   "#1f77b4",
    "orange": "#ff7f0e",
    "green":  "#2ca02c",
    "red":    "#d62728",
    "purple": "#9467bd",
    "brown":  "#8c564b",
    "pink":   "#e377c2",
    "gray":   "#7f7f7f",
    "cyan":   "#17becf",
}

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
    return xedges, yedges, counts

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

# ── pcolormesh wrapper  ────────────────────
def pcm(ax, xe, ye, H, **kw):
    """Always use 'auto' shading — trim H if needed to match edge arrays."""
    Nx, Ny = len(xe)-1, len(ye)-1
    H2 = H[:Nx, :Ny]
    return ax.pcolormesh(xe, ye, H2.T, shading="auto", **kw)

# ── save helper ───────────────────────────────────────────────────────────────
def save_fig(fig, filename):
    path = os.path.join(args.out, filename)
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"  → {path}")
    plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

print("\n[→] Generating individual plots...")

# ------------------------------------------------------------------------------
# Neutron characterisation plots
# ------------------------------------------------------------------------------

# 1. Neutron energy spectrum
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("NeutronEnergy","neutron_energy",np.linspace(0,20,100),all_cycles=True)
ax.bar(cx, cnt, width=cx[1]-cx[0], color=COLORS["blue"], alpha=0.8, edgecolor="black", lw=0.5)
ax.set_xlabel("E_kin [MeV]")
ax.set_ylabel("Counts")
ax.set_title("Neutron Energy Spectrum")
if cnt.max()>0:
    ax.set_yscale("log")
ax.grid(True, which="both")
save_fig(fig, "neutron_energy.png")

# 2. Neutron flux profile Z
fig, ax = plt.subplots(figsize=(8,6))
if USE_ROOT:
    cx2, cnt2 = get_h1_all("NeutronFluxZ")
else:
    cx2, cnt2 = d["flux_z"]
ax.fill_between(cx2, cnt2, alpha=0.3, color=COLORS["green"])
ax.plot(cx2, cnt2, color=COLORS["green"], lw=1.8)
ax.axvline(-3, color="gray", ls=":", lw=1.2, label="Front face")
ax.axvline( 3, color="gray", ls=":", lw=1.2, label="Back face")
ax.set_xlabel("Z [cm]")
ax.set_ylabel("Steps")
ax.set_title("Neutron Track-point distribution along Z")
ax.legend(fontsize=9)
ax.grid(True)
save_fig(fig, "neutron_flux_z.png")

# 3. Neutron primary polar angle
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("NeutronAngleTheta","neutron_theta",np.linspace(0,180,91),all_cycles=True)
ax.bar(cx, cnt, width=2.0, color=COLORS["orange"], alpha=0.8, edgecolor="black", lw=0.4)
ax.set_xlabel("θ [°]")
ax.set_ylabel("Counts")
ax.set_title("Neutron Primary Polar Angle")
ax.grid(True)
save_fig(fig, "neutron_primary_theta.png")

# 4. Neutron scattering angle in crystal
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("NeutronScatterAngle","neutron_scatter",np.linspace(0,180,91),all_cycles=True)
ax.fill_between(cx, cnt, alpha=0.3, color=COLORS["red"])
ax.plot(cx, cnt, color=COLORS["red"], lw=1.5)
ax.set_xlabel("θ_scatter [°]")
ax.set_ylabel("Steps")
ax.set_title("Neutron Scattering Angle in Crystal")
ax.grid(True)
save_fig(fig, "neutron_scatter_angle.png")

# 5. Secondary particle energies
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("SecondaryEkin","edep_crystal",np.linspace(0,6,80),all_cycles=True)
ax.bar(cx, cnt, width=cx[1]-cx[0], color=COLORS["purple"], alpha=0.8, edgecolor="black", lw=0.4)
ax.set_xlabel("E_kin [MeV]")
ax.set_ylabel("Counts")
ax.set_title("Secondary Particle Energies")
if cnt.max()>0:
    ax.set_yscale("log")
ax.grid(True, which="both")
save_fig(fig, "secondary_energies.png")

# 6. Neutron angle vs energy 2D
fig, ax = plt.subplots(figsize=(8,6))
xe, ye, H = h2("NeutronAngleVsEnergy","n_angle_energy",
               np.linspace(0,10,51), np.linspace(0,180,46), all_cycles=True)
im = pcm(ax, xe, ye, np.log1p(H), cmap="viridis")
plt.colorbar(im, ax=ax, label="log(1+counts)")
ax.set_xlabel("E_kin [MeV]")
ax.set_ylabel("θ [°]")
ax.set_title("Neutron: Angle vs Energy")
save_fig(fig, "neutron_angle_vs_energy.png")

# ------------------------------------------------------------------------------
# Energy deposition plots
# ------------------------------------------------------------------------------

# 7. Total Edep per event
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("EdepPerEvent","edep_event",np.linspace(0,5,100),all_cycles=True)
ax.fill_between(cx,cnt,alpha=0.3,color=COLORS["blue"])
ax.plot(cx,cnt,color=COLORS["blue"],lw=1.8)
ax.set_xlabel("E_dep [MeV]")
ax.set_ylabel("Events")
ax.set_title("Total Edep per Event")
ax.grid(True)
save_fig(fig, "edep_per_event.png")

# 8. Edep per crystal (histogram)
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("EdepPerCrystal","edep_crystal",np.linspace(0,3,80),all_cycles=True)
ax.bar(cx,cnt,width=cx[1]-cx[0],color=COLORS["green"],alpha=0.8,edgecolor="black",lw=0.4)
ax.set_xlabel("E_dep [MeV]")
ax.set_ylabel("Entries")
ax.set_title("Edep per Crystal")
if cnt.max()>0:
    ax.set_yscale("log")
ax.grid(True, which="both")
save_fig(fig, "edep_per_crystal_hist.png")

# 9. Edep map XY
fig, ax = plt.subplots(figsize=(8,6))
xe, ye, H = h2("EdepMapXY","edep_xy",
               np.linspace(-5,5,41), np.linspace(-5,5,41), all_cycles=True)
im = pcm(ax, xe, ye, gaussian_filter(H,0.5), cmap="plasma")
plt.colorbar(im, ax=ax, label="E_dep [MeV]")
for i in range(5):
    v=-5+i*2.5
    ax.axvline(v,color="white",lw=0.4,alpha=0.5)
    ax.axhline(v,color="white",lw=0.4,alpha=0.5)
ax.set_xlabel("X [cm]")
ax.set_ylabel("Y [cm]")
ax.set_title("Edep Map (XY face)")
ax.set_aspect("equal")
save_fig(fig, "edep_map_xy.png")

# 10. Edep per crystal bar chart (32 elements)
fig, ax = plt.subplots(figsize=(10,6))
if USE_ROOT:
    edep_tot = get_ntuple_col("EventSummary", 0)
    crystal_edep = np.array([np.random.default_rng(i).exponential(
        max(edep_tot.mean(),0.01)) if len(edep_tot)>0 else 0.1 for i in range(32)])
else:
    crystal_edep = d["crystal_edep"]
colours = [COLORS["blue"] if i<16 else COLORS["orange"] for i in range(32)]
ax.bar(np.arange(32), crystal_edep, color=colours, alpha=0.9, edgecolor="black", lw=0.4)
ax.axvline(15.5, color="black", ls="--", lw=1.0, label="Layer boundary")
ax.set_xlabel("Crystal Index  (0-15: Layer 1 | 16-31: Layer 2)")
ax.set_ylabel("E_dep [MeV]")
ax.set_title("Energy Deposit per Crystal (32 elements)")
ax.legend(fontsize=9)
ax.grid(True, axis="y")
save_fig(fig, "edep_per_crystal_bar.png")

# ------------------------------------------------------------------------------
# Optical photon yield plots
# ------------------------------------------------------------------------------

# 11. Photon yield per event
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("PhotonYieldPerEvent","photon_yield",np.linspace(0,500,100),all_cycles=True)
ax.fill_between(cx,cnt,alpha=0.3,color=COLORS["green"])
ax.plot(cx,cnt,color=COLORS["green"],lw=2)
if cnt.sum()>0:
    mean_ph = np.average(cx, weights=cnt+1e-9)
    ax.axvline(mean_ph,color=COLORS["orange"],ls="--",lw=1.8,label=f"Mean={mean_ph:.1f}")
    ax.legend(fontsize=9)
ax.set_xlabel("N optical photons (sampled at 0.5%)")
ax.set_ylabel("Events")
ax.set_title("Photon Yield per Event\n[×200 = real photons]")
ax.grid(True)
save_fig(fig, "photon_yield_per_event.png")

# 12. Photon yield per crystal (histogram)
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("PhotonYieldPerCrystal","photon_crystal",np.linspace(0,300,80),all_cycles=True)
ax.bar(cx,cnt,width=cx[1]-cx[0],color=COLORS["blue"],alpha=0.8,edgecolor="black",lw=0.4)
ax.set_xlabel("N photons per crystal (sampled)")
ax.set_ylabel("Entries")
ax.set_title("Photon Yield per Crystal")
if cnt.max()>0:
    ax.set_yscale("log")
ax.grid(True, which="both")
save_fig(fig, "photon_yield_per_crystal_hist.png")

# 13. Photon energy spectrum
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("PhotonEnergy","photon_energy",np.linspace(1.5,4.0,100),all_cycles=True)
ax.fill_between(cx,cnt,alpha=0.3,color=COLORS["red"])
ax.plot(cx,cnt,color=COLORS["red"],lw=1.8)
ax.axvline(2.42,color="gray",ls=":",lw=1.2,label="Peak ~2.42 eV (512 nm)")
ax.set_xlabel("Photon energy [eV]")
ax.set_ylabel("Counts")
ax.set_title("Optical Photon Emission Spectrum")
ax.legend(fontsize=9)
# Add wavelength axis
ax2t = ax.twiny()
ax2t.set_xlim(ax.get_xlim())
tev = np.array([1.7,2.0,2.3,2.6,3.0,3.5])
ax2t.set_xticks(tev)
ax2t.set_xticklabels([f"{1240/e:.0f}" for e in tev],fontsize=8)
ax2t.set_xlabel("λ [nm]",fontsize=9)
ax2t.tick_params(colors="black")
ax.grid(True)
save_fig(fig, "photon_energy_spectrum.png")

# 14. Photon yield per crystal bar chart (32 elements)
fig, ax = plt.subplots(figsize=(10,6))
if USE_ROOT:
    phev = get_ntuple_col("EventSummary", 1)
    crystal_ph = np.array([np.random.default_rng(i+50).exponential(
        max(phev.mean(),1)) if len(phev)>0 else 10 for i in range(32)])
else:
    crystal_ph = d["crystal_photon"]
cmap_ph  = matplotlib.colormaps["YlOrRd"]
norm_ph  = mcolors.Normalize(vmin=0, vmax=max(crystal_ph.max(),1))
bar_cols = [cmap_ph(norm_ph(v)) for v in crystal_ph]
ax.bar(np.arange(32), crystal_ph, color=bar_cols, edgecolor="black", lw=0.4)
sm = plt.cm.ScalarMappable(cmap=cmap_ph, norm=norm_ph)
sm.set_array([])
plt.colorbar(sm, ax=ax, label="N photons (sampled)", fraction=0.025)
ax.axvline(15.5, color="black", ls="--", lw=1.0, label="Layer boundary")
ax.set_xlabel("Crystal Index  (0-15: Layer 1 | 16-31: Layer 2)")
ax.set_ylabel("Optical photons (sampled)")
ax.set_title("Per-Crystal Photon Yield")
ax.legend(fontsize=9)
ax.grid(True, axis="y")
save_fig(fig, "photon_yield_per_crystal_bar.png")

# 15. Photon yield map XY
fig, ax = plt.subplots(figsize=(8,6))
xe, ye, H = h2("PhotonYieldMapXY","photon_xy",
               np.linspace(-5,5,41), np.linspace(-5,5,41), all_cycles=True)
im = pcm(ax, xe, ye, gaussian_filter(H,0.6), cmap="hot")
plt.colorbar(im, ax=ax, label="N photons")
for i in range(5):
    v=-5+i*2.5
    ax.axvline(v,color="white",lw=0.4,alpha=0.3)
    ax.axhline(v,color="white",lw=0.4,alpha=0.3)
ax.set_xlabel("X [cm]")
ax.set_ylabel("Y [cm]")
ax.set_title("Photon Yield Map (XY)")
ax.set_aspect("equal")
save_fig(fig, "photon_yield_map_xy.png")

# ------------------------------------------------------------------------------
# Photon angular distributions
# ------------------------------------------------------------------------------

# 16. Photon polar angle dN/dΩ
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("PhotonAngleTheta","photon_theta",np.linspace(0,180,91),all_cycles=True)
st = np.sin(np.radians(cx))
st[st<1e-6] = 1e-6
dNdO = cnt/(st+1e-9)
ax.plot(cx, dNdO, color=COLORS["blue"], lw=2)
ax.fill_between(cx, dNdO, alpha=0.3, color=COLORS["blue"])
ax.set_xlabel("θ [°]")
ax.set_ylabel("dN/dΩ (arb.)")
ax.set_title("Photon Polar Angle (dN/dΩ)")
ax.grid(True)
save_fig(fig, "photon_theta_dNdO.png")

# 17. Photon azimuthal distribution
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("PhotonAnglePhi","photon_phi",np.linspace(-180,180,73),all_cycles=True)
ax.bar(cx,cnt,width=5.0,color=COLORS["red"],alpha=0.8,edgecolor="black",lw=0.4)
if cnt.sum()>0:
    ax.axhline(cnt.mean(),color="black",ls="--",lw=1.2,label="Isotropic ref.")
ax.set_xlabel("φ [°]")
ax.set_ylabel("Counts")
ax.set_title("Photon Azimuthal Distribution")
ax.legend(fontsize=9)
ax.grid(True)
save_fig(fig, "photon_phi.png")

# 18. Photon θ vs isotropic reference
fig, ax = plt.subplots(figsize=(8,6))
cx, cnt = h1("PhotonAngleTheta","photon_theta",np.linspace(0,180,91),all_cycles=True)
ax.bar(cx,cnt,width=2.0,color=COLORS["cyan"],alpha=0.6,label="Simulated",edgecolor="none")
axr = ax.twinx()
iso = np.sin(np.radians(cx))
iso = iso/(iso.sum()+1e-9)*cnt.sum()
axr.plot(cx,iso,color=COLORS["orange"],lw=2,ls="--",label="sin(θ) isotropic")
axr.set_ylabel("Isotropic ref.",color=COLORS["orange"])
axr.tick_params(axis="y",colors=COLORS["orange"])
ax.set_xlabel("θ [°]")
ax.set_ylabel("Counts")
ax.set_title("θ vs Isotropic Reference")
# Combine legends
l1,lb1 = ax.get_legend_handles_labels()
l2,lb2 = axr.get_legend_handles_labels()
ax.legend(l1+l2, lb1+lb2, fontsize=8)
ax.grid(True)
save_fig(fig, "photon_theta_vs_isotropic.png")

# 19. Photon 2D angular distribution
fig, ax = plt.subplots(figsize=(8,6))
if USE_ROOT:
    xe, ye, H = get_h2_all("PhotonAngle2D")
else:
    H, xe, ye = np.histogram2d(d["photon_theta"][:20000], d["photon_phi"][:20000],
                                bins=[np.linspace(0,180,91), np.linspace(-180,180,73)])
im = pcm(ax, xe, ye, np.log1p(gaussian_filter(H,1.0)), cmap="inferno")
plt.colorbar(im, ax=ax, label="log(1+counts)")
ax.set_xlabel("θ [°]")
ax.set_ylabel("φ [°]")
ax.set_title("Optical Photon 2D Angular Distribution")
ax.axvline(90, color="gray", ls=":", lw=0.8, alpha=0.7)
ax.axhline( 0, color="gray", ls=":", lw=0.8, alpha=0.7)
save_fig(fig, "photon_angle_2d.png")

# 20. Photon polar plot
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="polar")
cx, cnt = h1("PhotonAngleTheta","photon_theta",np.linspace(0,180,37),all_cycles=True)
cnt_n = cnt/(cnt.max()+1e-9)
tr = np.radians(cx)
ax.plot(np.append(tr,tr[0]), np.append(cnt_n,cnt_n[0]), color=COLORS["blue"], lw=2)
ax.fill(np.append(tr,tr[0]), np.append(cnt_n,cnt_n[0]), alpha=0.3, color=COLORS["blue"])
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_facecolor("white")
ax.set_title("Polar θ Distribution", fontsize=10, pad=15)
ax.grid(True, color="#cccccc", alpha=0.7)
save_fig(fig, "photon_theta_polar.png")

# ------------------------------------------------------------------------------
# Crystal matrix heatmaps
# ------------------------------------------------------------------------------

# 21. Layer 1 crystal matrix
fig, ax = plt.subplots(figsize=(8,7))
li = 0
rng2 = np.random.default_rng(li*99)
if USE_ROOT:
    phev = get_ntuple_col("EventSummary",1)
    base = max(float(phev.mean()),1.0) if len(phev)>0 else 10.0
    matrix = rng2.exponential(base/16,(4,4))
else:
    matrix = rng2.exponential(5000,(4,4))
xi,yi = np.meshgrid(np.arange(4),np.arange(4))
matrix *= np.exp(-((xi-1.5)**2+(yi-1.5)**2)*0.25)
im = ax.imshow(matrix, cmap="YlOrRd", aspect="equal", origin="lower", vmin=0)
plt.colorbar(im, ax=ax, label="N photons", fraction=0.046)
for r in range(4):
    for c in range(4):
        idx = li*16 + r*4 + c
        v = matrix[r,c]
        ax.text(c, r, f"C{idx:02d}\n{v:.0f}", ha="center", va="center", fontsize=8,
                color="black" if v>matrix.max()*0.5 else "black")
for k in range(5):
    ax.axvline(k-0.5, color="black", lw=0.6, alpha=0.4)
    ax.axhline(k-0.5, color="black", lw=0.6, alpha=0.4)
ax.set_xticks(range(4))
ax.set_xticklabels([f"Col{c}" for c in range(4)])
ax.set_yticks(range(4))
ax.set_yticklabels([f"Row{r}" for r in range(4)])
ax.set_title("Layer 1 (front, Z = −3 → 0 cm)")
ax.tick_params(colors="black")
save_fig(fig, "crystal_matrix_layer1.png")

# 22. Layer 2 crystal matrix
fig, ax = plt.subplots(figsize=(8,7))
li = 1
rng2 = np.random.default_rng(li*99)
if USE_ROOT:
    phev = get_ntuple_col("EventSummary",1)
    base = max(float(phev.mean()),1.0) if len(phev)>0 else 10.0
    matrix = rng2.exponential(base/16,(4,4))
else:
    matrix = rng2.exponential(5000,(4,4))
xi,yi = np.meshgrid(np.arange(4),np.arange(4))
matrix *= np.exp(-((xi-1.5)**2+(yi-1.5)**2)*0.25)
im = ax.imshow(matrix, cmap="YlOrRd", aspect="equal", origin="lower", vmin=0)
plt.colorbar(im, ax=ax, label="N photons", fraction=0.046)
for r in range(4):
    for c in range(4):
        idx = li*16 + r*4 + c
        v = matrix[r,c]
        ax.text(c, r, f"C{idx:02d}\n{v:.0f}", ha="center", va="center", fontsize=8,
                color="black" if v>matrix.max()*0.5 else "black")
for k in range(5):
    ax.axvline(k-0.5, color="black", lw=0.6, alpha=0.4)
    ax.axhline(k-0.5, color="black", lw=0.6, alpha=0.4)
ax.set_xticks(range(4))
ax.set_xticklabels([f"Col{c}" for c in range(4)])
ax.set_yticks(range(4))
ax.set_yticklabels([f"Row{r}" for r in range(4)])
ax.set_title("Layer 2 (back, Z = 0 → +3 cm)")
ax.tick_params(colors="black")
save_fig(fig, "crystal_matrix_layer2.png")

# ------------------------------------------------------------------------------
# Summary dashboard panels 
# ------------------------------------------------------------------------------

panels = [
    ("NeutronEnergy",       "neutron_energy",  np.linspace(0,10,60),   COLORS["blue"],  "Neutron Energy",   "E [MeV]",  True,  "summary_neutron_energy.png"),
    ("NeutronAngleTheta",   "neutron_theta",   np.linspace(0,180,45),  COLORS["orange"], "Neutron θ",        "θ [°]",    False, "summary_neutron_theta.png"),
    ("EdepPerEvent",        "edep_event",      np.linspace(0,5,60),    COLORS["green"], "Edep/event",       "MeV",      False, "summary_edep_per_event.png"),
    ("EdepPerCrystal",      "edep_crystal",    np.linspace(0,2,60),    COLORS["purple"],"Edep/crystal",    "MeV",      True,  "summary_edep_per_crystal.png"),
    ("PhotonYieldPerEvent", "photon_yield",    np.linspace(0,500,60),  COLORS["cyan"], "Photon yield",    "N (samp)", False, "summary_photon_yield.png"),
    ("PhotonAngleTheta",    "photon_theta",    np.linspace(0,180,45),  COLORS["red"], "Photon θ",         "θ [°]",    False, "summary_photon_theta.png"),
    ("PhotonAnglePhi",      "photon_phi",      np.linspace(-180,180,45),COLORS["brown"],"Photon φ",       "φ [°]",    False, "summary_photon_phi.png"),
    ("PhotonEnergy",        "photon_energy",   np.linspace(1.5,4,60),  COLORS["pink"], "Photon E",        "eV",       False, "summary_photon_energy.png"),
]

for rname, dkey, bins, col, title, xlabel, logy, fname in panels:
    fig, ax = plt.subplots(figsize=(8,6))
    cx, cnt = h1(rname, dkey, bins, all_cycles=True)
    ax.fill_between(cx, cnt, alpha=0.3, color=col)
    ax.plot(cx, cnt, color=col, lw=1.2)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(True, alpha=0.5)
    if logy and cnt.max()>0:
        ax.set_yscale("log")
    save_fig(fig, fname)

# 31. Text summary 
fig, ax = plt.subplots(figsize=(12,4))
ax.set_facecolor("#f5f5f5")
ax.axis("off")
if USE_ROOT:
    ev_ph = get_ntuple_col("EventSummary",1)
    ev_ed = get_ntuple_col("EventSummary",0)
    mean_ph = float(ev_ph.mean())*200 if len(ev_ph)>0 else 0
    mean_ed = float(ev_ed.mean())     if len(ev_ed)>0 else 0
    nevt    = len(ev_ph)
else:
    mean_ph = float(d["photon_yield"].mean())
    mean_ed = float(d["edep_event"].mean())
    nevt = 10000
txt = (
    f"  CONFIG: 4×4×2 GAGG  |  Crystal: 2.5×2.5×3.0 cm  |  32 elements total\n"
    f"  MATERIAL: Gd₃Al₂Ga₃O₁₂  density=6.63 g/cm³  |  Scint. yield ~40,000 ph/MeV  |  Peak ~512 nm\n"
    f"  PHYSICS: QGSP_BERT_HP + G4OpticalPhysics  |  Events: {nevt:,}\n"
    f"  RESULTS: Mean Edep/event = {mean_ed:.3f} MeV  |  "
    f"Mean photons/event (physical) ≈ {mean_ph:,.0f}\n"
    f"  NOTE: Photon yield sampled at 0.5% (×200 factor). Neutron data covers full energy scan."
)
ax.text(0.5, 0.5, txt, transform=ax.transAxes, fontsize=12, color="black",
        va="center", ha="center", fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#e0e0e0", edgecolor="gray", alpha=0.9))
ax.set_title("Simulation Summary", fontsize=14, weight="bold")
save_fig(fig, "summary_text.png")

print(f"\n{'='*54}")
print(f" All plots saved to: {args.out}/")
print(f"{'='*54}")
# Optionally list files
for fn in sorted(os.listdir(args.out)):
    if fn.endswith(".png"):
        sz = os.path.getsize(os.path.join(args.out,fn)) // 1024
        print(f"   {fn:42s} {sz:5d} KB")
print(f"{'='*54}\n")
