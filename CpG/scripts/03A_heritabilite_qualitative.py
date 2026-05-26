"""
Pipeline HPC - Etape 3A : Heritabilite qualitative (M/U binaire)
Methode : Correlation de Pearson entre qual_P0 et qual_R2
Seuils testes : 5 et 10 reads
"""
import os, time
import polars as pl
import numpy as np
from scipy import stats
from tqdm import tqdm
import pandas as pd

IN_DIR  = os.path.expanduser("~/Documents/temporaire/results/filtrage")
OUT_DIR = os.path.expanduser("~/Documents/temporaire/results/heritabilite")
os.makedirs(OUT_DIR, exist_ok=True)

SEUILS = [5, 10]

INDIVIDUS = [
    "IT_01_01","IT_01_03","IT_01_05","IT_01_06",
    "IT_02_01","IT_02_03","IT_02_05","IT_02_07",
    "IT_03_01","IT_03_03","IT_03_05","IT_03_06",
    "IT_04_01","IT_04_03","IT_04_05","IT_04_06",
    "IT_05_01","IT_05_03","IT_05_05","IT_05_06",
    "IT_06_01","IT_06_02","IT_06_04","IT_06_06",
    "IT_07_02","IT_07_04","IT_07_05","IT_07_07",
]

# Populations pour grouper les resultats
POPULATIONS = {
    "IT_01": ["IT_01_01","IT_01_03","IT_01_05","IT_01_06"],
    "IT_02": ["IT_02_01","IT_02_03","IT_02_05","IT_02_07"],
    "IT_03": ["IT_03_01","IT_03_03","IT_03_05","IT_03_06"],
    "IT_04": ["IT_04_01","IT_04_03","IT_04_05","IT_04_06"],
    "IT_05": ["IT_05_01","IT_05_03","IT_05_05","IT_05_06"],
    "IT_06": ["IT_06_01","IT_06_02","IT_06_04","IT_06_06"],
    "IT_07": ["IT_07_02","IT_07_04","IT_07_05","IT_07_07"],
}

def calculer_h2_qualitatif(nom, seuil):
    """
    Calcule h² qualitatif par correlation de Pearson
    entre qual_P0 (binaire) et qual_R2 (binaire).
    h² = r² (coefficient de determination)
    """
    filepath = os.path.join(IN_DIR, f"{nom}_seuil{seuil}.parquet")
    df = pl.read_parquet(filepath).select(["site","chr","qual_P0","qual_R2"])

    # Conversion en numpy pour scipy
    q_p0 = df["qual_P0"].to_numpy()
    q_r2 = df["qual_R2"].to_numpy()

    # Correlation de Pearson globale
    r, pval = stats.pearsonr(q_p0, q_r2)
    h2_global = round(r**2, 6)

    # Proportion de sites methyles
    pct_M_P0 = round(q_p0.mean() * 100, 2)
    pct_M_R2 = round(q_r2.mean() * 100, 2)

    # Tableau de contingence pour verifier
    # MM / MU / UM / UU
    MM = int(((q_p0 == 1) & (q_r2 == 1)).sum())
    MU = int(((q_p0 == 1) & (q_r2 == 0)).sum())
    UM = int(((q_p0 == 0) & (q_r2 == 1)).sum())
    UU = int(((q_p0 == 0) & (q_r2 == 0)).sum())

    # h² par chromosome
    h2_par_chr = []
    for chr_name, group in df.group_by("chr"):
        gp0 = group["qual_P0"].to_numpy()
        gr2 = group["qual_R2"].to_numpy()
        if len(gp0) < 100:
            continue
        if gp0.std() == 0 or gr2.std() == 0:
            continue
        r_chr, p_chr = stats.pearsonr(gp0, gr2)
        h2_par_chr.append({
            "individu"  : nom,
            "chr"       : chr_name[0],
            "seuil"     : seuil,
            "n_sites"   : len(gp0),
            "h2"        : round(r_chr**2, 6),
            "r"         : round(r_chr, 6),
            "pval"      : round(p_chr, 6),
        })

    return {
        "global": {
            "individu"  : nom,
            "seuil"     : seuil,
            "n_sites"   : len(df),
            "h2"        : h2_global,
            "r"         : round(r, 6),
            "pval"      : round(pval, 10),
            "pct_M_P0"  : pct_M_P0,
            "pct_M_R2"  : pct_M_R2,
            "MM"        : MM,
            "MU"        : MU,
            "UM"        : UM,
            "UU"        : UU,
        },
        "par_chr": h2_par_chr
    }

if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"  Etape 3A : Heritabilite QUALITATIVE (M/U binaire)")
    print(f"  Methode  : Correlation Pearson — h² = r²")
    print(f"  Seuils   : {SEUILS} reads")
    print(f"{'='*60}\n")

    t_total = time.time()
    tous_global  = []
    tous_par_chr = []

    for seuil in SEUILS:
        print(f"\n--- Seuil {seuil} reads ---")
        for nom in tqdm(INDIVIDUS, desc=f"h² seuil{seuil}",
                        unit="individu", ncols=65):
            res = calculer_h2_qualitatif(nom, seuil)
            tous_global.append(res["global"])
            tous_par_chr.extend(res["par_chr"])
            tqdm.write(
                f"  {nom} | h²={res['global']['h2']:.4f} "
                f"| %M_P0={res['global']['pct_M_P0']}% "
                f"| %M_R2={res['global']['pct_M_R2']}%"
            )

    # Sauvegarder
    df_global  = pd.DataFrame(tous_global).sort_values(["seuil","individu"])
    df_par_chr = pd.DataFrame(tous_par_chr).sort_values(["seuil","individu","chr"])

    df_global.to_csv(os.path.join(OUT_DIR, "h2_qualitatif_global.csv"), index=False)
    df_par_chr.to_csv(os.path.join(OUT_DIR, "h2_qualitatif_par_chr.csv"), index=False)

    duree_totale = round(time.time() - t_total, 1)

    # Résumé par seuil
    print(f"\n{'='*60}")
    for seuil in SEUILS:
        sub = df_global[df_global["seuil"] == seuil]
        print(f"\nSeuil {seuil} reads :")
        print(f"  h² moyen    : {sub['h2'].mean():.4f}")
        print(f"  h² min      : {sub['h2'].min():.4f}")
        print(f"  h² max      : {sub['h2'].max():.4f}")
        print(f"  %M P0 moyen : {sub['pct_M_P0'].mean():.1f}%")
        print(f"  %M R2 moyen : {sub['pct_M_R2'].mean():.1f}%")

    print(f"\nTemps total : {duree_totale}s ({round(duree_totale/60,1)} min)")
    print(f"Resultats   : {OUT_DIR}")
    print(f"{'='*60}")
