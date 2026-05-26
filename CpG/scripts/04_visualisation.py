"""
Pipeline HPC - Etape 4 : Visualisation des resultats
Graphiques : h² qual vs quant, heatmap chromosomes,
             scatter P0 vs R2, comparaison seuils
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import polars as pl

HER_DIR  = os.path.expanduser("~/Documents/temporaire/results/heritabilite")
FILT_DIR = os.path.expanduser("~/Documents/temporaire/results/filtrage")
OUT_DIR  = os.path.expanduser("~/Documents/temporaire/results/figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Style global
plt.rcParams.update({
    'font.family'   : 'DejaVu Sans',
    'font.size'     : 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi'    : 150,
})

POPULATIONS = {
    "IT_01": ["IT_01_01","IT_01_03","IT_01_05","IT_01_06"],
    "IT_02": ["IT_02_01","IT_02_03","IT_02_05","IT_02_07"],
    "IT_03": ["IT_03_01","IT_03_03","IT_03_05","IT_03_06"],
    "IT_04": ["IT_04_01","IT_04_03","IT_04_05","IT_04_06"],
    "IT_05": ["IT_05_01","IT_05_03","IT_05_05","IT_05_06"],
    "IT_06": ["IT_06_01","IT_06_02","IT_06_04","IT_06_06"],
    "IT_07": ["IT_07_02","IT_07_04","IT_07_05","IT_07_07"],
}
POP_COLORS = {
    "IT_01":"#4C72B0","IT_02":"#DD8452","IT_03":"#55A868",
    "IT_04":"#C44E52","IT_05":"#8172B3","IT_06":"#937860",
    "IT_07":"#DA8BC3"
}

def get_pop(ind):
    for pop, membres in POPULATIONS.items():
        if ind in membres:
            return pop
    return "unknown"

# Charger les donnees
qual_g  = pd.read_csv(os.path.join(HER_DIR, "h2_qualitatif_global.csv"))
quant_g = pd.read_csv(os.path.join(HER_DIR, "h2_quantitatif_global.csv"))
qual_c  = pd.read_csv(os.path.join(HER_DIR, "h2_qualitatif_par_chr.csv"))
quant_c = pd.read_csv(os.path.join(HER_DIR, "h2_quantitatif_par_chr.csv"))

qual_g["population"]  = qual_g["individu"].apply(get_pop)
quant_g["population"] = quant_g["individu"].apply(get_pop)

print("Creation des figures...")

# ─────────────────────────────────────────────────
# FIGURE 1 : h² qual vs quant par individu (seuil 5)
# ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle("Heritabilite CpG — Qualitatif vs Quantitatif (seuil 5 reads)",
             fontsize=14, fontweight='bold', y=1.01)

for ax, seuil in zip(axes, [5, 10]):
    q5  = qual_g[qual_g["seuil"] == seuil].sort_values("individu")
    qt5 = quant_g[quant_g["seuil"] == seuil].sort_values("individu")

    x     = np.arange(len(q5))
    width = 0.35
    colors_ind = [POP_COLORS[get_pop(i)] for i in q5["individu"]]

    bars1 = ax.bar(x - width/2, q5["h2"],  width, label="Qualitatif (M/U)",
                   color=colors_ind, alpha=0.7, edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, qt5["h2"], width, label="Quantitatif (%)",
                   color=colors_ind, alpha=1.0, edgecolor='black', linewidth=0.5)

    ax.set_xlabel("Individu")
    ax.set_ylabel("h²")
    ax.set_title(f"Seuil {seuil} reads")
    ax.set_xticks(x)
    ax.set_xticklabels(q5["individu"], rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0.88, 1.0)
    ax.axhline(y=q5["h2"].mean(),  color='steelblue', linestyle='--',
               alpha=0.7, linewidth=1, label=f"Moy qual={q5['h2'].mean():.3f}")
    ax.axhline(y=qt5["h2"].mean(), color='darkorange', linestyle='--',
               alpha=0.7, linewidth=1, label=f"Moy quant={qt5['h2'].mean():.3f}")

    # Legende populations
    pop_patches = [mpatches.Patch(color=c, label=p)
                   for p, c in POP_COLORS.items()]
    leg1 = ax.legend(handles=pop_patches, title="Population",
                     loc='upper left', fontsize=7, ncol=2)
    ax.add_artist(leg1)
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig1_h2_qual_vs_quant.png"),
            bbox_inches='tight')
plt.close()
print("  fig1 OK")

# ─────────────────────────────────────────────────
# FIGURE 2 : Heatmap h² par chromosome (quant seuil 5)
# ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.suptitle("Heritabilite quantitative par chromosome",
             fontsize=14, fontweight='bold')

for ax, seuil in zip(axes, [5, 10]):
    sub = quant_c[quant_c["seuil"] == seuil].copy()
    # Garder chromosomes principaux
    chrs_principaux = [f"Fvb{i}" for i in range(1, 8)]
    sub = sub[sub["chr"].isin(chrs_principaux)]
    pivot = sub.pivot_table(index="individu", columns="chr",
                            values="h2", aggfunc="mean")
    pivot = pivot[chrs_principaux]

    sns.heatmap(pivot, ax=ax, cmap="YlOrRd", vmin=0.93, vmax=1.0,
                annot=True, fmt=".3f", linewidths=0.5,
                cbar_kws={"label": "h²"}, annot_kws={"size": 7})
    ax.set_title(f"h² quantitatif — seuil {seuil} reads")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Individu")
    ax.tick_params(axis='x', rotation=0)
    ax.tick_params(axis='y', rotation=0, labelsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig2_heatmap_chromosomes.png"),
            bbox_inches='tight')
plt.close()
print("  fig2 OK")

# ─────────────────────────────────────────────────
# FIGURE 3 : Scatter P0 vs R2 pour 4 individus
# ─────────────────────────────────────────────────
exemples = ["IT_01_01","IT_02_01","IT_05_05","IT_06_02"]
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("Methylation P0 vs R2 — Approche quantitative (seuil 5)",
             fontsize=14, fontweight='bold')

for ax, nom in zip(axes.flatten(), exemples):
    df = pl.read_parquet(
        os.path.join(FILT_DIR, f"{nom}_seuil5.parquet")
    ).select(["methyl_P0","methyl_R2"]).sample(n=50000, seed=42).to_pandas()

    h2_val = quant_g[(quant_g["individu"]==nom) &
                     (quant_g["seuil"]==5)]["h2"].values[0]
    sl_val = quant_g[(quant_g["individu"]==nom) &
                     (quant_g["seuil"]==5)]["slope"].values[0]
    pop    = get_pop(nom)

    ax.scatter(df["methyl_P0"], df["methyl_R2"],
               alpha=0.05, s=1, color=POP_COLORS[pop])

    # Ligne de regression
    x_line = np.linspace(0, 100, 100)
    ic_val = quant_g[(quant_g["individu"]==nom) &
                     (quant_g["seuil"]==5)]["intercept"].values[0]
    ax.plot(x_line, sl_val * x_line + ic_val,
            color='red', linewidth=1.5, label=f"y={sl_val:.3f}x+{ic_val:.2f}")
    ax.plot([0,100],[0,100], 'k--', linewidth=1, alpha=0.4, label="y=x")

    ax.set_xlabel("% Methylation P0 (parent)")
    ax.set_ylabel("% Methylation R2 (clone gen.3)")
    ax.set_title(f"{nom} ({pop})\nh²={h2_val:.4f} | slope={sl_val:.4f}")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig3_scatter_P0_vs_R2.png"),
            bbox_inches='tight')
plt.close()
print("  fig3 OK")

# ─────────────────────────────────────────────────
# FIGURE 4 : Comparaison h² par population + seuil
# ─────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Comparaison h² par population et seuil de couverture",
             fontsize=14, fontweight='bold')

datasets = [
    (qual_g,  "Qualitatif",  "steelblue"),
    (quant_g, "Quantitatif", "darkorange"),
]

for row, (df, label, color) in enumerate(datasets):
    for col, seuil in enumerate([5, 10]):
        ax = axes[row][col]
        sub = df[df["seuil"] == seuil].copy()
        sub["population"] = sub["individu"].apply(get_pop)

        pop_h2 = sub.groupby("population")["h2"].agg(["mean","std"]).reset_index()
        pop_h2 = pop_h2.sort_values("population")

        bars = ax.bar(pop_h2["population"], pop_h2["mean"],
                      yerr=pop_h2["std"], capsize=4,
                      color=[POP_COLORS[p] for p in pop_h2["population"]],
                      edgecolor='black', linewidth=0.5, alpha=0.85)

        ax.set_xlabel("Population")
        ax.set_ylabel("h² moyen")
        ax.set_title(f"{label} — seuil {seuil} reads\n"
                     f"h² global = {sub['h2'].mean():.4f} ± {sub['h2'].std():.4f}")
        ax.set_ylim(0.88, 1.0)
        ax.axhline(y=sub["h2"].mean(), color='black',
                   linestyle='--', linewidth=1, alpha=0.6)
        ax.grid(axis='y', alpha=0.3)

        # Valeurs sur les barres
        for bar, val in zip(bars, pop_h2["mean"]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                    f"{val:.3f}", ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4_h2_par_population.png"),
            bbox_inches='tight')
plt.close()
print("  fig4 OK")

# ─────────────────────────────────────────────────
# FIGURE 5 : Boxplot h² qual vs quant tous seuils
# ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
fig.suptitle("Distribution h² : Qualitatif vs Quantitatif",
             fontsize=14, fontweight='bold')

data_box = []
for seuil in [5, 10]:
    for _, row in qual_g[qual_g["seuil"]==seuil].iterrows():
        data_box.append({"Methode":"Qualitatif",
                         "Seuil":f"seuil {seuil}", "h2":row["h2"]})
    for _, row in quant_g[quant_g["seuil"]==seuil].iterrows():
        data_box.append({"Methode":"Quantitatif",
                         "Seuil":f"seuil {seuil}", "h2":row["h2"]})

df_box = pd.DataFrame(data_box)
sns.boxplot(data=df_box, x="Seuil", y="h2", hue="Methode",
            palette={"Qualitatif":"steelblue","Quantitatif":"darkorange"},
            ax=ax, width=0.5)
sns.stripplot(data=df_box, x="Seuil", y="h2", hue="Methode",
              palette={"Qualitatif":"steelblue","Quantitatif":"darkorange"},
              ax=ax, dodge=True, size=4, alpha=0.6, legend=False)

ax.set_ylabel("h²")
ax.set_xlabel("")
ax.set_ylim(0.88, 1.0)
ax.grid(axis='y', alpha=0.3)
ax.legend(title="Methode", fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig5_boxplot_comparaison.png"),
            bbox_inches='tight')
plt.close()
print("  fig5 OK")

print(f"\nToutes les figures sauvegardees dans : {OUT_DIR}")
print("Fichiers :")
for f in sorted(os.listdir(OUT_DIR)):
    print(f"  {f}")
