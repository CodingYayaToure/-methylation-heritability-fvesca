#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script : 05_delta_par_population.py
But    : Visualiser le gain de l'approche quantitative (delta = h²_quant - h²_qual)
         en regroupant les 28 individus par population (IT_01 à IT_07).
Sorties : - delta_par_population_barplot.png (moyenne + barre d'erreur)
         - delta_par_population_boxplot.png (boîte à moustaches)
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Chemins (à adapter si nécessaire)
HER_DIR = os.path.expanduser("~/Documents/temporaire/results/heritabilite")
OUT_DIR = os.path.expanduser("~/Documents/temporaire/results/figures")
os.makedirs(OUT_DIR, exist_ok=True)

# 1. Charger les données qualitatives et quantitatives (seuil 5 reads)
qual = pd.read_csv(os.path.join(HER_DIR, "h2_qualitatif_global.csv"))
quant = pd.read_csv(os.path.join(HER_DIR, "h2_quantitatif_global.csv"))

# Ne conserver que le seuil 5 (ou changer si vous préférez seuil 10)
qual5 = qual[qual["seuil"] == 5][["individu", "h2"]].rename(columns={"h2": "h2_qual"})
quant5 = quant[quant["seuil"] == 5][["individu", "h2"]].rename(columns={"h2": "h2_quant"})

# 2. Fusionner et calculer le delta
merged = pd.merge(qual5, quant5, on="individu")
merged["delta"] = merged["h2_quant"] - merged["h2_qual"]

# 3. Ajouter la population (IT_01, IT_02, ...)
merged["population"] = merged["individu"].str[:5]  # "IT_01" à "IT_07"

# Vérification rapide
print(merged.head())

# 4. Statistiques par population
pop_stats = merged.groupby("population")["delta"].agg(["mean", "std", "count"]).reset_index()
pop_stats = pop_stats.sort_values("population")
print(pop_stats)

# 5. Graphique 1 : barres avec barres d'erreur (écart-type)
plt.figure(figsize=(8, 5))
plt.bar(pop_stats["population"], pop_stats["mean"], yerr=pop_stats["std"],
        capsize=5, color="steelblue", edgecolor="black", alpha=0.8)
plt.ylabel("Gain de l'approche quantitative\n(h²_quant - h²_qual)")
plt.xlabel("Population")
plt.title("Héritabilité CpG (seuil 5 reads) : avantage de l'approche quantitative")
plt.axhline(y=0, color="red", linestyle="--", linewidth=1)
plt.grid(axis="y", alpha=0.3)
# Ajout des valeurs au-dessus des barres
for i, row in pop_stats.iterrows():
    plt.text(i, row["mean"] + 0.002, f"{row['mean']:.3f}", ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "delta_par_population_barplot.png"), dpi=150)
plt.close()
print("Figure sauvegardée : delta_par_population_barplot.png")

# 6. Graphique 2 : boxplot (plus détaillé)
plt.figure(figsize=(8, 5))
sns.boxplot(data=merged, x="population", y="delta", palette="Set2", showfliers=False)
sns.stripplot(data=merged, x="population", y="delta", color="black", alpha=0.5, size=4)
plt.ylabel("Delta h² (quant - qual)")
plt.xlabel("Population")
plt.title("Distribution du gain apporté par l'approche quantitative")
plt.axhline(y=0, color="red", linestyle="--", linewidth=1)
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "delta_par_population_boxplot.png"), dpi=150)
plt.close()
print("Figure sauvegardée : delta_par_population_boxplot.png")

print(f"Toutes les figures sont dans : {OUT_DIR}")
