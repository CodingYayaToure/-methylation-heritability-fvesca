"""
Pipeline HPC - Etape 3B : Heritabilite quantitative (% methylation)
Methode : Correlation de Pearson sur valeurs continues 0-100%
         + Regression lineaire pour visualisation
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

POPULATIONS = {
    "IT_01": ["IT_01_01","IT_01_03","IT_01_05","IT_01_06"],
    "IT_02": ["IT_02_01","IT_02_03","IT_02_05","IT_02_07"],
    "IT_03": ["IT_03_01","IT_03_03","IT_03_05","IT_03_06"],
    "IT_04": ["IT_04_01","IT_04_03","IT_04_05","IT_04_06"],
    "IT_05": ["IT_05_01","IT_05_03","IT_05_05","IT_05_06"],
    "IT_06": ["IT_06_01","IT_06_02","IT_06_04","IT_06_06"],
    "IT_07": ["IT_07_02","IT_07_04","IT_07_05","IT_07_07"],
}

def calculer_h2_quantitatif(nom, seuil):
    """
    Calcule h² quantitatif par correlation de Pearson
    entre methyl_P0 (%) et methyl_R2 (%).
    h² = r² (coefficient de determination)
    Aussi : regression lineaire slope + intercept
    """
    filepath = os.path.join(IN_DIR, f"{nom}_seuil{seuil}.parquet")
    df = pl.read_parquet(filepath).select(
        ["site","chr","methyl_P0","methyl_R2"]
    )

    m_p0 = df["methyl_P0"].to_numpy()
    m_r2 = df["methyl_R2"].to_numpy()

    # Correlation de Pearson
    r, pval = stats.pearsonr(m_p0, m_r2)
    h2 = round(r**2, 6)

    # Regression lineaire y = ax + b
    slope, intercept, r_reg, p_reg, se = stats.linregress(m_p0, m_r2)

    # Statistiques descriptives
    mean_P0 = round(float(np.mean(m_p0)), 3)
    mean_R2 = round(float(np.mean(m_r2)), 3)
    std_P0  = round(float(np.std(m_p0)),  3)
    std_R2  = round(float(np.std(m_r2)),  3)

    # h² par chromosome
    h2_par_chr = []
    for chr_name, group in df.group_by("chr"):
        gp0 = group["methyl_P0"].to_numpy()
        gr2 = group["methyl_R2"].to_numpy()
        if len(gp0) < 100:
            continue
        if gp0.std() == 0 or gr2.std() == 0:
            continue
        r_chr, p_chr = stats.pearsonr(gp0, gr2)
        sl, ic, _, _, _ = stats.linregress(gp0, gr2)
        h2_par_chr.append({
            "individu" : nom,
            "chr"      : chr_name[0],
            "seuil"    : seuil,
            "n_sites"  : len(gp0),
            "h2"       : round(r_chr**2, 6),
            "r"        : round(r_chr, 6),
            "slope"    : round(sl, 6),
            "intercept": round(ic, 6),
            "pval"     : round(p_chr, 10),
        })

    return {
        "global": {
            "individu"  : nom,
            "seuil"     : seuil,
            "n_sites"   : len(df),
            "h2"        : h2,
            "r"         : round(r, 6),
            "pval"      : round(pval, 10),
            "slope"     : round(slope, 6),
            "intercept" : round(intercept, 6),
            "mean_P0"   : mean_P0,
            "mean_R2"   : mean_R2,
            "std_P0"    : std_P0,
            "std_R2"    : std_R2,
        },
        "par_chr": h2_par_chr
    }

if __name__ == "__main__":
    print(f"{'='*60}")
    print(f"  Etape 3B : Heritabilite QUANTITATIVE (% methylation)")
    print(f"  Methode  : Correlation Pearson — h² = r²")
    print(f"             + Regression lineaire slope/intercept")
    print(f"  Seuils   : {SEUILS} reads")
    print(f"{'='*60}\n")

    t_total = time.time()
    tous_global  = []
    tous_par_chr = []

    for seuil in SEUILS:
        print(f"\n--- Seuil {seuil} reads ---")
        for nom in tqdm(INDIVIDUS, desc=f"h² seuil{seuil}",
                        unit="individu", ncols=65):
            res = calculer_h2_quantitatif(nom, seuil)
            tous_global.append(res["global"])
            tous_par_chr.extend(res["par_chr"])
            tqdm.write(
                f"  {nom} | h²={res['global']['h2']:.4f} "
                f"| slope={res['global']['slope']:.4f} "
                f"| mean_P0={res['global']['mean_P0']:.1f}% "
                f"| mean_R2={res['global']['mean_R2']:.1f}%"
            )

    # Sauvegarder
    df_global  = pd.DataFrame(tous_global).sort_values(["seuil","individu"])
    df_par_chr = pd.DataFrame(tous_par_chr).sort_values(["seuil","individu","chr"])

    df_global.to_csv(
        os.path.join(OUT_DIR, "h2_quantitatif_global.csv"), index=False)
    df_par_chr.to_csv(
        os.path.join(OUT_DIR, "h2_quantitatif_par_chr.csv"), index=False)

    duree_totale = round(time.time() - t_total, 1)

    # Résumé par seuil
    print(f"\n{'='*60}")
    for seuil in SEUILS:
        sub = df_global[df_global["seuil"] == seuil]
        print(f"\nSeuil {seuil} reads :")
        print(f"  h² moyen    : {sub['h2'].mean():.4f}")
        print(f"  h² min      : {sub['h2'].min():.4f}")
        print(f"  h² max      : {sub['h2'].max():.4f}")
        print(f"  slope moyen : {sub['slope'].mean():.4f}")
        print(f"  mean_P0 moy : {sub['mean_P0'].mean():.1f}%")
        print(f"  mean_R2 moy : {sub['mean_R2'].mean():.1f}%")

    print(f"\nTemps total : {duree_totale}s ({round(duree_totale/60,1)} min)")
    print(f"Resultats   : {OUT_DIR}")
    print(f"{'='*60}")
