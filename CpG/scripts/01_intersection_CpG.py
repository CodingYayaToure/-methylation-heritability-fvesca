"""
Pipeline HPC - Intersection P0/R2 CpG
Moteur : Polars (parallélisme natif multi-thread)
"""
import os, time, gzip
import polars as pl
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

DATA_DIR = os.path.expanduser("~/Documents/temporaire/CpG")
OUT_DIR  = os.path.expanduser("~/Documents/temporaire/results")
os.makedirs(OUT_DIR, exist_ok=True)

INDIVIDUS = [
    "FV_IT_01_01","FV_IT_01_03","FV_IT_01_05","FV_IT_01_06",
    "FV_IT_02_01","FV_IT_02_03","FV_IT_02_05","FV_IT_02_07",
    "FV_IT_03_01","FV_IT_03_03","FV_IT_03_05","FV_IT_03_06",
    "FV_IT_04_01","FV_IT_04_03","FV_IT_04_05","FV_IT_04_06",
    "FV_IT_05_01","FV_IT_05_03","FV_IT_05_05","FV_IT_05_06",
    "FV_IT_06_01","FV_IT_06_02","FV_IT_06_04","FV_IT_06_06",
    "FV_IT_07_02","FV_IT_07_04","FV_IT_07_05","FV_IT_07_07",
]

PAIRES = [
    (f"{ind}_P0_WC0_M1_CpG.bedGraph.gz",
     f"{ind}_R2_GG0_M1_CpG.bedGraph.gz",
     ind.replace("FV_",""))
    for ind in INDIVIDUS
]

def lire_bedgraph_polars(filepath):
    """Lecture ultra-rapide avec Polars — parallélisme natif."""
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
        .with_columns(
            (pl.col("chr") + "_" + pl.col("start").cast(pl.Utf8)).alias("site")
        )
    )

def traiter_paire(args):
    """Intersection P0 ∩ R2 via Polars join (multi-thread natif)."""
    fichier_P0, fichier_R2, nom = args
    t0 = time.time()

    p0 = lire_bedgraph_polars(os.path.join(DATA_DIR, fichier_P0))
    r2 = lire_bedgraph_polars(os.path.join(DATA_DIR, fichier_R2))

    # Join inner = intersection rapide
    merged = (
        p0.select(["site","chr","start","methyl_pct"])
          .rename({"methyl_pct": "methyl_P0"})
          .join(
              r2.select(["site","methyl_pct"])
                .rename({"methyl_pct": "methyl_R2"}),
              on="site",
              how="inner"
          )
          .with_columns([
              (pl.col("methyl_P0") >= 50).cast(pl.Int8).alias("qual_P0"),
              (pl.col("methyl_R2") >= 50).cast(pl.Int8).alias("qual_R2"),
          ])
    )

    # Sauvegarde parquet
    out = os.path.join(OUT_DIR, f"{nom}_intersection.parquet")
    merged.write_parquet(out)

    duree = round(time.time() - t0, 1)
    return {
        "individu"      : nom,
        "sites_P0"      : len(p0),
        "sites_R2"      : len(r2),
        "sites_communs" : len(merged),
        "pct_commun"    : round(len(merged) / min(len(p0), len(r2)) * 100, 2),
        "duree_s"       : duree,
    }

if __name__ == "__main__":
    import multiprocessing
    n_cpu = multiprocessing.cpu_count()

    print(f"{'='*55}")
    print(f"  Pipeline HPC — Polars + ProcessPoolExecutor")
    print(f"  CPUs : {n_cpu} | Paires : {len(PAIRES)}")
    print(f"{'='*55}\n")

    t_total = time.time()
    resultats = []

    # ProcessPoolExecutor pour lancer les paires en parallèle
    # Polars utilise déjà tous les threads en interne par paire
    n_workers = min(n_cpu, 4)
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(traiter_paire, p): p[2] for p in PAIRES}
        with tqdm(total=len(PAIRES), desc="Progression",
                  unit="paire", ncols=65) as pbar:
            for future in as_completed(futures):
                res = future.result()
                resultats.append(res)
                pbar.set_postfix({"dernier": res["individu"],
                                  "sites": f"{res['sites_communs']:,}"})
                pbar.update(1)

    # Résumé
    import pandas as pd
    resume = pd.DataFrame(resultats).sort_values("individu")
    duree_totale = round(time.time() - t_total, 1)

    print(f"\n{'='*55}")
    print(resume.to_string(index=False))
    print(f"\nTemps total  : {duree_totale}s ({round(duree_totale/60,1)} min)")
    print(f"Résultats    : {OUT_DIR}")
    print(f"{'='*55}")

    resume.to_csv(os.path.join(OUT_DIR, "resume_intersection.csv"), index=False)
