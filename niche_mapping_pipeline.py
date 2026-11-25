"""
Tumour-immune niche mapping in cutaneous squamous cell carcinoma
Dataset: GEO GSE144240 (Ji et al., Cell 2020)
  GSE144239 - spatial transcriptomics (12 legacy ST sections, 4 patients)
  GSE144236 - matched scRNA-seq reference (48,164 cells)

Pipeline
  1. Fetch ST counts (spot coordinates are encoded in the column names) and the scRNA reference.
  2. Build an 18-type cell-type signature by streaming the scRNA matrix gene-by-gene
     (avoids holding 48k cells in memory).
  3. Deconvolve each spot by Poisson EM (the RCTD likelihood, without per-gene overdispersion).
  4. Cluster spots into niches on own + neighbourhood composition, Harmony-integrated across patients.
  5. Permutation tests for niche adjacency and ligand-receptor signalling at the interface.

Requires: numpy, pandas, scipy, scanpy, harmonypy==0.0.9 (newer releases pull in torch)
"""
import urllib.request, gzip, os
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc, anndata as ad
from scipy.stats import spearmanr

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/"
os.makedirs("data", exist_ok=True)

# ---------------------------------------------------------------- 1. download
for url, path in [
    (GEO + "GSE144239/suppl/GSE144239_ST_all_counts.txt.gz", "data/st_counts.txt.gz"),
    (GEO + "GSE144236/suppl/GSE144236_cSCC_counts.txt.gz",   "data/sc_counts.txt.gz"),
    (GEO + "GSE144236/suppl/GSE144236_patient_metadata_new.txt.gz", "data/sc_meta.txt.gz")]:
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)

# ST matrix: rows = Patient / Section replicate / pixel_x / pixel_y, then genes.
# Column names carry the array coordinates, e.g. "P2_7x45_1" = patient P2, x=7, y=45, replicate 1.
st = pd.read_csv("data/st_counts.txt.gz", sep="\t", index_col=0)
obs = pd.DataFrame({
    "patient": "P" + st.loc["Patient"].astype(int).astype(str).values,
    "rep":     st.loc["Section replicate"].astype(int).values,
    "pixel_x": st.loc["pixel_x"].astype(float).values,
    "pixel_y": st.loc["pixel_y"].astype(float).values}, index=st.columns.str.strip('"'))
obs["section"] = obs.patient + "_rep" + obs.rep.astype(str)
xy = pd.Series(obs.index).str.extract(r"^P\d+_(\d+)x(\d+)_\d+$").astype(float)
obs["array_x"], obs["array_y"] = xy[0].values, xy[1].values
genes = list(st.index[4:])
X = sp.csr_matrix(st.iloc[4:].to_numpy(dtype=np.float32).T)      # spots x genes
del st

# ------------------------------------------------- 2. cell-type signature matrix
# 24 level-2 labels collapsed to 18 deconvolution targets: the four tumour
# keratinocyte states stay separate, normal keratinocytes merge, immune stays resolved.
COLLAPSE = {"Tumor_KC_Basal":"Tumor_KC_Basal","Tumor_KC_Cyc":"Tumor_KC_Cyc","Tumor_KC_Diff":"Tumor_KC_Diff",
    "TSK":"TSK","Keratinocyte":"Normal_KC","Normal_KC_Basal":"Normal_KC","Normal_KC_Diff":"Normal_KC",
    "Normal_KC_Cyc":"Normal_KC","Pilosebaceous":"Skin_appendage","Eccrine":"Skin_appendage",
    "Melanocyte":"Melanocyte","Fibroblast":"Fibroblast","Endothelial Cell":"Endothelial","Tcell":"T_cell",
    "NK":"NK","B Cell":"B_cell","Mac":"Macrophage","MDSC":"MDSC","LC":"Langerhans","CD1C":"cDC2_CD1C",
    "CLEC9A":"cDC1_CLEC9A","ASDC":"pDC_ASDC","PDC":"pDC_ASDC"}

md = pd.read_csv("data/sc_meta.txt.gz", sep="\t", index_col=0)
with gzip.open("data/sc_counts.txt.gz", "rt") as fh:
    cells = [c.strip('"') for c in fh.readline().rstrip("\n").split("\t")]
m = md.reindex(cells)
lab = m.level2_celltype.map(COLLAPSE)
keep = lab.notna().values & (m["tum.norm"] == "Tumor").values     # tumour tissue, multiplets dropped
cat = pd.Categorical(lab[keep]); types = list(cat.categories); codes = cat.codes.astype(np.int64)
scale = (1e4 / m.loc[keep, "nCount_RNA"].values)                  # per-cell CPM-10k scaling
kidx = np.where(keep)[0]

rows, sums = [], []
with gzip.open("data/sc_counts.txt.gz", "rt") as fh:
    for _ in range(3): fh.readline()                              # header + 2 metadata rows
    for line in fh:
        i = line.index("\t"); g = line[:i].strip('"')
        v = np.fromstring(line[i+1:], sep="\t")[kidx] * scale
        rows.append(g); sums.append(np.bincount(codes, weights=v, minlength=len(types)))
sig = pd.DataFrame(np.vstack(sums) / np.bincount(codes, minlength=len(types)), index=rows, columns=types)
sig = sig.groupby(level=0).mean()

# ------------------------------------------------------- 3. Poisson EM deconvolution
# y_g ~ Poisson(sum_k n_k * s_gk); s_gk = mean counts of gene g in one cell of type k.
# Multiplicative EM update keeps n_k >= 0 and, unlike Gaussian NNLS on these 18
# collinear keratinocyte profiles, does not collapse onto a single winner.
common = [g for g in genes if g in sig.index]
gidx = {g: i for i, g in enumerate(genes)}
S = sig.loc[common]
spec = S.div(S.sum(1) + 1e-9, axis=0)                             # per-gene specificity
sel = sorted({g for ct in S.columns for g in spec[ct][S[ct] > 0.5].sort_values(ascending=False).index[:250]})
Sm = S.loc[sel].to_numpy(np.float64) + 1e-6
Y  = np.asarray(X[:, [gidx[g] for g in sel]].todense(), dtype=np.float64)
n  = np.ones((Y.shape[0], Sm.shape[1])); colsum = Sm.sum(0)
for _ in range(300):
    n *= ((Y / (n @ Sm.T + 1e-9)) @ Sm) / colsum
prop = pd.DataFrame(n / n.sum(1, keepdims=True), index=obs.index, columns=S.columns)

# validation: deconvolved fraction vs canonical marker expression in the same spot
tot = np.asarray(X.sum(1)).ravel()
cpm = np.asarray(X.todense(), dtype=np.float32) / (tot[:, None] + 1e-9) * 1e4
mscore = lambda gs: cpm[:, [gidx[g] for g in gs if g in gidx]].sum(1)
for ct, gs in {"T_cell": ["CD3D","CD3E","CD2","TRBC2"], "Macrophage": ["CD68","CD163","C1QA","C1QB"],
               "Fibroblast": ["COL1A1","COL1A2","DCN","LUM"], "Endothelial": ["PECAM1","VWF","CLDN5"],
               "B_cell": ["MS4A1","CD79A","IGHM"]}.items():
    print(f"  {ct}: rho = {spearmanr(prop[ct], mscore(gs)).statistic:.2f}")

# -------------------------------------------------------------- 4. niche clustering
obs["n_counts"] = tot; obs["n_genes"] = np.asarray((X > 0).sum(1)).ravel()

def neighbour_graph(obs, section):
    """8-connected neighbours on the ST array grid (spacing 1 unit)."""
    idx = obs.index[obs.section == section]
    xy = obs.loc[idx, ["array_x", "array_y"]].values
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
    return obs.index.get_indexer(idx), (d <= 1.5) & (d > 0)

smooth = np.zeros_like(prop.values)
for s in obs.section.unique():
    ii, A = neighbour_graph(obs, s)
    Aself = A + np.eye(len(ii), dtype=bool)
    smooth[ii] = (Aself @ prop.values[ii]) / Aself.sum(1, keepdims=True)

# features: own composition + neighbourhood composition, sqrt for variance stabilisation
F = np.hstack([np.sqrt(prop.values), np.sqrt(smooth)]).astype(np.float32)
a = ad.AnnData(F, obs=obs.copy())
sc.pp.scale(a); sc.tl.pca(a, n_comps=20)
# without integration the niches are patient-specific (mixing entropy ~0.2 of a max of 1)
sc.external.pp.harmony_integrate(a, "patient", basis="X_pca", adjusted_basis="X_h")
sc.pp.neighbors(a, n_neighbors=15, use_rep="X_h")
sc.tl.leiden(a, resolution=0.4, key_added="niche", flavor="igraph", n_iterations=2, directed=False)
obs["niche"] = a.obs.niche.values
comp = prop.groupby(obs.niche.values, observed=True).mean()

# ------------------------------------------- 5a. niche adjacency permutation test
labs = sorted(obs.niche.astype(str).unique()); K = len(labs)
code = pd.Categorical(obs.niche.astype(str), categories=labs).codes
rng = np.random.default_rng(1); NPERM = 500
pairs = np.zeros((K, K)); null = np.zeros((NPERM, K, K))
for s in obs.section.unique():
    ii, A = neighbour_graph(obs, s)
    r, c = np.nonzero(np.triu(A)); cs = code[ii]
    np.add.at(pairs, (cs[r], cs[c]), 1); np.add.at(pairs, (cs[c], cs[r]), 1)
    for p in range(NPERM):                       # labels permuted within section
        pc = rng.permutation(cs)
        np.add.at(null[p], (pc[r], pc[c]), 1); np.add.at(null[p], (pc[c], pc[r]), 1)
adj_z = pd.DataFrame((pairs - null.mean(0)) / (null.std(0) + 1e-9), index=labs, columns=labs)

# ------------------------------- 5b. cell-type co-localisation (within-patient)
coloc = sum(spearmanr(prop.loc[idx].values).statistic
            for idx in obs.groupby("patient", observed=True).groups.values()) / obs.patient.nunique()
coloc = pd.DataFrame(coloc, index=prop.columns, columns=prop.columns)

# ------------------------------- 5c. ligand-receptor scoring across the interface
# Curated TME pair list - the ST platform (~1.6k counts/spot) is too shallow for a
# full database screen, so this is a targeted panel, not an unbiased one.
LR = [("CXCL9","CXCR3"),("CXCL10","CXCR3"),("CXCL11","CXCR3"),("CXCL12","CXCR4"),("CCL19","CCR7"),
      ("CCL21","CCR7"),("CCL2","CCR2"),("CCL5","CCR5"),("CCL17","CCR4"),("CCL22","CCR4"),("CXCL1","CXCR2"),
      ("CXCL8","CXCR1"),("CSF1","CSF1R"),("IL34","CSF1R"),("CD274","PDCD1"),("CD86","CTLA4"),
      ("LGALS9","HAVCR2"),("TGFB1","TGFBR1"),("TGFB1","TGFBR2"),("IL10","IL10RA"),("IL6","IL6R"),
      ("VEGFA","FLT1"),("SPP1","CD44"),("THBS1","CD47"),("MIF","CD74"),("ANXA1","FPR1"),
      ("LAMC2","ITGB4"),("TNC","ITGB1"),("IL1B","IL1R1"),("PVR","TIGIT"),("CD70","CD27")]

E = ad.AnnData(X.astype(np.float32), obs=obs.copy(), var=pd.DataFrame(index=genes))
E.var_names_make_unique(); sc.pp.normalize_total(E, target_sum=1e4); sc.pp.log1p(E)
M = E.X.toarray(); gi = {g: i for i, g in enumerate(E.var_names)}
LR = [(l, r) for l, r in LR if l in gi and r in gi]

def interface_z(src_niches, dst_niches, nperm=300, seed=2):
    """Ligand in src spot x receptor in adjacent dst spot, vs randomly drawn adjacent pairs."""
    nic = obs.niche.astype(str).values
    pr, pc = [], []
    for s in obs.section.unique():
        ii, A = neighbour_graph(obs, s)
        r, c = np.nonzero(A); pr.append(ii[r]); pc.append(ii[c])
    pr, pc = np.concatenate(pr), np.concatenate(pc)
    mask = np.isin(nic[pr], src_niches) & np.isin(nic[pc], dst_niches)
    L = np.array([gi[l] for l, _ in LR]); R = np.array([gi[r] for _, r in LR])
    score = (M[pr[mask]][:, L] * M[pc[mask]][:, R]).mean(0)
    rng = np.random.default_rng(seed)
    nulls = np.vstack([(M[pr[p]][:, L] * M[pc[p]][:, R]).mean(0)
                       for p in (rng.permutation(len(pr))[:mask.sum()] for _ in range(nperm))])
    return pd.DataFrame({"ligand": [l for l, _ in LR], "receptor": [r for _, r in LR],
                         "score": score, "z": (score - nulls.mean(0)) / (nulls.std(0) + 1e-9),
                         "n_pairs": mask.sum()})

# niche ids below are the Leiden labels from this run - re-derive them from `comp`
# (TSK-high/immune-low vs fibroblast+immune-high) before reusing on new data.
pd.concat([obs, prop], axis=1).to_csv("spot_niche_table.csv")
comp.to_csv("niche_cell_composition.csv")
adj_z.to_csv("niche_adjacency_zscores.csv")
coloc.to_csv("celltype_colocalisation.csv")
