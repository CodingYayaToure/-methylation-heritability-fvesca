"""
Pipeline HPC - Etape 2 : Filtrage par couverture minimum
Traitement sequentiel pour eviter saturation RAM
Polars gere le parallelisme en interne par fichier
"""
import os, time
import polars as pl
from tqdm import tqdm
import pandas as pd

DATA_DIR_CpG = os.path.expanduser("~/Documents/temporaire/CpG")
IN_DIR       = os.path.expanduser("~/Documents/temporaire/results")
OUT_DIR      = os.path.expanduser("~/Documents/temporaire/results/filtrage")
os.makedirs(OUT_DIR, exist_ok=True)

SEUILS = [5, 10, 15, 20]

INDIVIDUS = [
    "IT_01_01","IT_01_03","IT_01_05","IT_01_06",
    "IT_02_01","IT_02_03","IT_02_05","IT_02_07",
    "IT_03_01","IT_03_03","IT_03_05","IT_03_06",
    "IT_04_01","IT_04_03","IT_04_05","IT_04_06",
    "IT_05_01","IT_05_03","IT_05_05","IT_05_06",
    "IT_06_01","IT_06_02","IT_06_04","IT_06_06",
    "IT_07_02","IT_07_04","IT_07_05","IT_07_07",
]

def lire_coverage(filepath):
    """Lit uniquement site + coverage pour economiser la RAM."""
    return (
        pl.read_csv(
            filepath,
            separator="\t",
            has_header=False,
            comment_prefix="t",
            new_columns=["chr","start","end","methyl_pct","reads_M","reads_U"],
            schema_overrides={
                "chr"       : pl.Utf8,
                "start"     : pl.Int32,
                "end"       : pl.Int32,
                "methyl_pct": pl.Float32,
                "reads_M"   : pl.Int32,
                "reads_U"   : pl.Int32,
            },
            ignore_errors=True,
        )
        .with_columns([
            (pl.col("chr") + "_" + pl.col("start").cast(pl.Utf8)).alias("site"),
            (pl.col("reads_M") + pl.col("reads_U")).alias("coverage"),
        ])
        .select(["site", "coverage"])
    )

def filtrer_individu(nom):
    t0 = time.time()

    # Charger intersection
    inter = pl.read_parquet(os.path.join(IN_DIR, f"{nom}_intersection.parquet"))

    # Charger seulement site+coverage de P0 et R2
    fv_id  = "FV_" + nom
    cov_p0 = lire_coverage(
        os.path.join(DATA_DIR_CpG, f"{fv_id}_P0_WC0_M1_CpG.bedGraph.gz"))
    cov_r2 = lire_coverage(
        os.path.join(DATA_DIR_CpG, f"{fv_id}_R2_GG0_M1_CpG.bedGraph.gz"))

    # Joindre les coverages
    merged = (
        inter
        .join(cov_p0.rename({"coverage":"cov_P0"}), on="site", how="left")
        .join(cov_r2.rename({"coverage":"cov_R2"}), on="site", how="left")
        .fill_null(0)
    )

    resultats_seuils = []
    for seuil in SEUILS:
        filtre = (
            merged
            .filter(
                (pl.col("cov_P0") >= seuil) &
                (pl.col("cov_R2") >= seuil)
            )
            .with_columns([
                (pl.col("methyl_P0") >= 50).cast(pl.Int8).alias("qual_P0"),
                (pl.col("methyl_R2") >= 50).cast(pl.Int8).alias("qual_R2"),
            ])
        )

        out = os.path.join(OUT_DIR, f"{nom}_seuil{seuil}.parquet")
        filtre.write_parquet(out)

        resultats_seuils.append({
            "individu"     : nom,
            "seuil"        : seuil,
            "sites_avant"  : len(merged),
            "sites_apres"  : len(filtre),
            "pct_conserve" : round(len(filtre) / len(merged) * 100, 2),
        })

    duree = round(time.time() - t0, 1)
    for r in resultats_seuils:
        r["duree_s"] = duree

    # Liberer la memoire explicitement
    del inter, cov_p0, cov_r2, merged

    return resultats_seuils

if __name__ == "__main__":
    import multiprocessing
    n_cpu = multiprocessing.cpu_count()

    print(f"{'='*60}")
    print(f"  Pipeline HPC — Etape 2 : Filtrage couverture")
    print(f"  Seuils : {SEUILS} reads | Individus : {len(INDIVIDUS)}")
    print(f"  Mode   : sequentiel (Polars multi-thread interne)")
    print(f"{'='*60}\n")

    t_total = time.time()
    tous_resultats = []

    for nom in tqdm(INDIVIDUS, desc="Filtrage", unit="individu", ncols=65):
        res = filtrer_individu(nom)
        tous_resultats.extend(res)
        tqdm.write(f"  {nom} -> seuil5: {res[0]['sites_apres']:,} sites ({res[0]['pct_conserve']}%)")

    resume = pd.DataFrame(tous_resultats).sort_values(["seuil","individu"])
    duree_totale = round(time.time() - t_total, 1)

    print(f"\n{'='*60}")
    for seuil in SEUILS:
        sub = resume[resume["seuil"] == seuil]
        moy_sites  = int(sub["sites_apres"].mean())
        moy_conserve = sub["pct_conserve"].mean()
        print(f"Seuil {seuil:2d} reads | "
              f"Sites moyens : {moy_sites:>9,} | "
              f"Conserves : {moy_conserve:.1f}%")

    print(f"\nTemps total : {duree_totale}s ({round(duree_totale/60,1)} min)")
    print(f"Resultats   : {OUT_DIR}")
    print(f"{'='*60}")

    resume.to_csv(os.path.join(OUT_DIR, "resume_filtrage.csv"), index=False)
