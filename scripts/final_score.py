"""Final scoring + top-10 CIF collection. Runs inside Docker (has atomworks)."""
import json, glob, shutil
from pathlib import Path
import numpy as np

PIPELINE   = Path("/workspace")
BOLTZ_BASE = PIPELINE / "boltz_outputs/boltz_results_boltz_inputs/predictions"
RFD3_OUT   = PIPELINE / "outputs/rfd3"
MPNN_OUT   = PIPELINE / "outputs/mpnn"
FINAL      = PIPELINE / "final"
STRUCTS    = PIPELINE / "top_structures"
FINAL.mkdir(exist_ok=True); STRUCTS.mkdir(exist_ok=True)

from atomworks.io.utils.io_utils import load_any
import biotite.structure.io.pdbx as pdbx

def kabsch_rmsd(P, Q):
    P=P-P.mean(0); Q=Q-Q.mean(0)
    U,S,Vt=np.linalg.svd(P.T@Q); d=np.linalg.det(Vt.T@U.T)
    R=Vt.T@np.diag([1,1,d])@U.T
    return float(np.sqrt(np.mean(np.sum((P@R.T-Q)**2,axis=1))))

def boltz_ca(cif,chain):
    aa=pdbx.get_structure(pdbx.CIFFile.read(str(cif)),model=1,include_bonds=False)
    m=(aa.atom_name=="CA")&(aa.chain_id==chain); return aa.coord[m]

def rfd3_ca(cif,r1,r2):
    raw=load_any(str(cif)); aa=raw[0] if hasattr(raw,"__getitem__") else raw
    ch=list(set(aa.chain_id))[0]
    m=(aa.chain_id==ch)&np.isin(aa.atom_name,["CA"])&(aa.res_id>=r1)&(aa.res_id<=r2)
    return aa.coord[m]

mpnn_data=json.load(open(MPNN_OUT/"all_sequences.json"))
mpnn_data.sort(key=lambda x:x["sequence_recovery"])
seen_bb,top50=[],[]
for r in mpnn_data:
    if r["backbone"] not in seen_bb:
        seen_bb.append(r["backbone"]); top50.append(r)
    if len(top50)>=50: break

pred_dirs={d.name:d for d in BOLTZ_BASE.iterdir() if d.is_dir()}

results=[]
for idx,r in enumerate(top50):
    bb,mb,rec=r["backbone"],r["minibinder_sequence"],r["sequence_recovery"]
    rank=idx+1; name=f"rank{rank:02d}_{bb}_rec{rec:.3f}"
    pred_dir=pred_dirs.get(name)
    plddt=iptm_ab=iptm_bc=0.0; sc_rmsd=999.0
    if pred_dir:
        conf=pred_dir/f"confidence_{name}_model_0.json"
        cif =pred_dir/f"{name}_model_0.cif"
        if conf.exists():
            c=json.load(open(conf))
            plddt=float(c.get("complex_plddt",0))
            p=c.get("pair_chains_iptm",{})
            iptm_ab=float(p.get("0",{}).get("1",0) or p.get("1",{}).get("0",0))
            iptm_bc=float(p.get("1",{}).get("2",0) or p.get("2",{}).get("1",0))
        if cif.exists():
            try:
                rc=RFD3_OUT/f"{bb}.cif"
                if rc.exists():
                    rv=rfd3_ca(rc,1,115); rm=rfd3_ca(rc,116,170); rn=rfd3_ca(rc,171,285)
                    bv=boltz_ca(cif,"A"); bm=boltz_ca(cif,"B"); bn=boltz_ca(cif,"C")
                    nv=min(len(bv),len(rv)); nn=min(len(bn),len(rn)); nm=min(len(bm),len(rm))
                    if nv>10 and nn>10 and nm>5:
                        P=np.vstack([bv[:nv],bn[:nn]]); Q=np.vstack([rv[:nv],rn[:nn]])
                        Pc=P-P.mean(0); Qc=Q-Q.mean(0)
                        U,S,Vt=np.linalg.svd(Pc.T@Qc); d=np.linalg.det(Vt.T@U.T)
                        Rmat=Vt.T@np.diag([1,1,d])@U.T
                        sc_rmsd=kabsch_rmsd((bm[:nm]-P.mean(0))@Rmat.T+Q.mean(0),rm[:nm])
            except Exception as e: print(f"  scRMSD err {bb}: {e}")
    passed=sc_rmsd<2.0 and plddt>0.60 and iptm_ab>0.3 and iptm_bc>0.3
    results.append({"rank":rank,"backbone":bb,"minibinder_sequence":mb,
                    "sequence_recovery":rec,"boltz_plddt":round(plddt,3),
                    "boltz_plddt_pct":round(plddt*100,1),
                    "iptm_MB_VH1":round(iptm_ab,3),"iptm_MB_Nb":round(iptm_bc,3),
                    "sc_rmsd_A":round(sc_rmsd,3),"pass_filter":passed})

results.sort(key=lambda x:(x["sc_rmsd_A"] if x["sc_rmsd_A"]<900 else 999,
                             -(x["iptm_MB_VH1"]*x["iptm_MB_Nb"])**0.5))
json.dump(results,open(FINAL/"final_results.json","w"),indent=2)

# Top-10 CIF collection
scored=[(r,(r["iptm_MB_VH1"]*r["iptm_MB_Nb"])**0.5) for r in results if r["iptm_MB_VH1"]>0]
scored.sort(key=lambda x:-x[1])
top10=[r for r,_ in scored[:10]]

print("\n=== Top-10 by sqrt(ipTM_VH1 x ipTM_Nb) ===")
for i,r in enumerate(top10):
    bb=r["backbone"]; rec=r["sequence_recovery"]; rank=r["rank"]
    combined=(r["iptm_MB_VH1"]*r["iptm_MB_Nb"])**0.5
    name=f"rank{rank:02d}_{bb}_rec{rec:.3f}"
    pred_dir=pred_dirs.get(name)

    rfd3_src=RFD3_OUT/f"{bb}.cif"
    if rfd3_src.exists():
        shutil.copy(rfd3_src, STRUCTS/f"top{i+1:02d}_{bb}_rfd3.cif")

    if pred_dir:
        boltz_src=pred_dir/f"{name}_model_0.cif"
        if boltz_src.exists():
            shutil.copy(boltz_src, STRUCTS/f"top{i+1:02d}_{bb}_boltz2.cif")

    print(f"  top{i+1:02d} {bb}  sqrt_ipTM={combined:.3f}  pLDDT={r['boltz_plddt_pct']:.0f}%  "
          f"ipTM_VH1={r['iptm_MB_VH1']:.3f}  ipTM_Nb={r['iptm_MB_Nb']:.3f}  scRMSD={r['sc_rmsd_A']:.1f}A")

with open(STRUCTS/"top10_summary.tsv","w") as f:
    f.write("rank\tbackbone\tcombined_iptm\tiptm_VH1\tiptm_Nb\tpLDDT_pct\tsc_rmsd_A\tminibinder_sequence\n")
    for i,r in enumerate(top10):
        c=(r["iptm_MB_VH1"]*r["iptm_MB_Nb"])**0.5
        f.write(f"{i+1}\t{r['backbone']}\t{c:.3f}\t{r['iptm_MB_VH1']:.3f}\t"
                f"{r['iptm_MB_Nb']:.3f}\t{r['boltz_plddt_pct']:.1f}\t{r['sc_rmsd_A']:.2f}\t"
                f"{r['minibinder_sequence']}\n")

print(f"\nSaved {len(list(STRUCTS.glob('*.cif')))} CIF files + summary TSV")
print(f"\n{'='*76}")
print(f"{'R':>3}  {'Backbone':>28}  {'scRMSD':>7}  {'pLDDT%':>6}  {'VH1':>5}  {'Nb':>5}")
print("  "+"-"*74)
for r in results[:20]:
    print(f"{r['rank']:>3}.  {r['backbone']:>28}  {r['sc_rmsd_A']:>7.3f}  "
          f"{r['boltz_plddt_pct']:>6.1f}  {r['iptm_MB_VH1']:>5.3f}  {r['iptm_MB_Nb']:>5.3f}")
print(f"\nPassed: {sum(1 for r in results if r['pass_filter'])}/{len(results)}")
