# Analyse de l'héritabilité — Contexte CpG

Ce dossier contient l'ensemble du pipeline et des résultats pour l'étude de l'héritabilité
de la méthylation des dinucléotides CpG chez *Fragaria vesca* (fraisier des bois).

*Partie du projet [methylation-heritability-fvesca](../README.md) · Stage M1 CHPS UPVD 2025–2026*

---

## Seuil de binarisation (approche qualitative)

Pour le contexte CpG, le seuil utilisé dans l'approche qualitative (binaire) est **50 %**.
Ce seuil correspond à la moyenne du pourcentage de méthylation observée dans les données (~52,9 %),
arrondie à 50 %. Il est défini dans `scripts/03A_heritabilite_qualitative.py` via `THRESHOLD = 50`.

---

## Organisation du dossier CpG/

```
CpG/
├── README.md
├── data/                        <- fichiers bruts .bedGraph.gz (non versionnés, 56 fichiers)
├── scripts/
│   ├── 01_intersection_CpG.py
│   ├── 02_filtrage_couverture.py
│   ├── 03A_heritabilite_qualitative.py
│   ├── 03B_heritabilite_quantitative.py
│   ├── 04_visualisation.py
│   └── 05_comparaison_delta.py
└── results/
    ├── figures/
    │   ├── fig1_h2_qual_vs_quant.png
    │   ├── fig2_heatmap_chromosomes.png
    │   ├── fig3_scatter_P0_vs_R2.png
    │   ├── fig4_h2_par_population.png
    │   ├── fig5_boxplot_comparaison.png
    │   ├── delta_par_population_barplot.png
    │   ├── delta_par_population_boxplot.png
    │   └── delta_par_individu.png
    ├── heritabilite/
    │   ├── h2_qualitatif_global.csv
    │   ├── h2_qualitatif_par_chr.csv
    │   ├── h2_quantitatif_global.csv
    │   └── h2_quantitatif_par_chr.csv
    ├── filtrage/
    │   └── resume_filtrage.csv
    ├── resume_intersection.csv
    └── individual_information.csv
```

---

## Description des fichiers d'entrée (data/)

Les fichiers bruts sont au format **bedGraph** compressé (`.gz`).
Chaque fichier correspond à un individu et une génération clonale :

- `P0` : génération parentale
- `R2` : clone de 3e génération (fille)

Nom typique : `FV_IT_01_01_P0_WC0_M1_CpG.bedGraph.gz` (population `IT_01`, individu `01`, génération `P0`)

| Colonne | Nom | Type | Description |
|---------|-----|------|-------------|
| 1 | `chrom` | chaîne | Chromosome (Fvb1 à Fvb7) |
| 2 | `start` | entier | Position de début (1-based) |
| 3 | `end` | entier | Position de fin = start + 1 |
| 4 | `meth_pct` | entier | Pourcentage de méthylation (0–100) |
| 5 | `reads_m` | entier | Nombre de reads méthylés |
| 6 | `reads_um` | entier | Nombre de reads non méthylés |

Couverture = `reads_m + reads_um`. Seuil minimal utilisé : **5 reads**.

---

## Pipeline d'analyse (6 étapes)

| Étape | Script | Description | Sortie principale |
|-------|--------|-------------|-------------------|
| 1 | `01_intersection_CpG.py` | Sites communs P0 ∩ R2 par individu | `results/IT_*_intersection.parquet` |
| 2 | `02_filtrage_couverture.py` | Filtrage couverture (5, 10, 15, 20 reads) | `results/filtrage/*_seuilX.parquet` |
| 3A | `03A_heritabilite_qualitative.py` | Binarisation à 50 % puis h² (Pearson r²) | `results/heritabilite/h2_qualitatif_*.csv` |
| 3B | `03B_heritabilite_quantitative.py` | h² sur valeurs continues + régression | `results/heritabilite/h2_quantitatif_*.csv` |
| 4 | `04_visualisation.py` | Génération de toutes les figures | `results/figures/` |
| 5 | `05_comparaison_delta.py` | Delta h² par individu et par population | `results/figures/delta_*.png` |

Commande unique depuis `CpG/` :

```bash
make all
```

Pour une étape isolée : `make intersection`, `make filtrage`, etc. (voir le `Makefile`).

---

## Résultats numériques (seuil de couverture 5 reads)

| Approche | h² moyen (± écart-type) | Interprétation |
|----------|------------------------|----------------|
| **Quantitative** | **0,9714 ± 0,004** | Héritabilité quasi-totale, très stable entre individus |
| Qualitative | 0,9342 ± 0,017 | Sous-estime de 3,7 %, 4 fois plus variable |

Gain apporté par l'approche quantitative : de +0,032 (IT_01) à +0,047 (IT_07).
Les populations de haute altitude (IT_06, IT_07) présentent un delta plus élevé,
suggérant plus de sites CpG avec des pourcentages intermédiaires proches du seuil.

---

## Figures

### Figure 1 — h² qualitatif vs quantitatif par individu

![fig1](results/figures/fig1_h2_qual_vs_quant.png)

**Statistique** : barres pour chaque individu (28) — bleu = qualitatif, orange = quantitatif —
aux seuils de couverture 5 et 10. La barre orange est systématiquement plus haute.

**Biologie** : La méthylation CpG est transmise avec une fidélité exceptionnelle d'une génération
clonale à l'autre. L'approche quantitative révèle une héritabilité encore plus forte car elle tient
compte de l'intensité de la méthylation (ex. 51 % vs 90 %), alors que la binarisation les
place dans la même catégorie.

---

### Figure 2 — Heatmap de l'héritabilité par chromosome

![fig2](results/figures/fig2_heatmap_chromosomes.png)

**Statistique** : valeurs de h² quantitatif par chromosome (Fvb1 à Fvb7) pour chaque individu ;
toutes > 0,95, avec une homogénéité remarquable.

**Biologie** : La stabilité de la méthylation n'est pas confinée à quelques régions ;
elle est généralisable à l'ensemble du génome. Les mécanismes de maintenance
(méthyltransférases) agissent uniformément.

---

### Figure 3 — Nuage de points P0 vs R2

![fig3](results/figures/fig3_scatter_P0_vs_R2.png)

**Statistique** : 50 000 sites aléatoires pour 4 individus représentatifs.
Pente de la régression ≈ 1, intercept ≈ 0, r² ≈ 0,97.

**Biologie** : Chaque point est un site CpG : sa valeur chez le parent (P0) est quasi-identique
chez le clone de 3e génération (R2). L'épigénome est clonable avec une très haute fidélité.
Les écarts correspondent à de rares épimutations spontanées ou à du bruit technique.

---

### Figure 4 — h² moyen par population

![fig4](results/figures/fig4_h2_par_population.png)

**Statistique** : h² moyen par population (IT_01 à IT_07) pour qualitatif et quantitatif.
Les barres quantitatives sont toujours supérieures, aucune population ne montre de baisse.

**Biologie** : L'héritabilité épigénétique est indépendante de l'altitude d'origine.
Les clones de haute altitude conservent leur paysage méthylé avec la même précision.
C'est un trait espèce-spécifique, non altéré par l'environnement passé.

---

### Figure 5 — Distribution comparée (boxplot)

![fig5](results/figures/fig5_boxplot_comparaison.png)

**Statistique** : distribution des h² (qual vs quant). La médiane du quantitatif est plus haute
et sa dispersion est beaucoup plus faible. L'écart-type passe de 0,017 à 0,004.

**Biologie** : L'approche quantitative donne une estimation plus précise et reproductible.
La forte variance du qualitatif reflète l'arbitraire du seuil : deux individus avec des
distributions légèrement différentes autour de 50 % donneront des h² qualitatifs très
différents, alors que leur héritabilité réelle est presque identique.

---

### Figure 6 — Delta h² par population

![delta barplot](results/figures/delta_par_population_barplot.png)

![delta boxplot](results/figures/delta_par_population_boxplot.png)

**Statistique** : gain `h²_quant − h²_qual` moyenné par population.
IT_07 : 0,047 ; IT_03 : 0,044 ; IT_01 : 0,032.
Variabilité intra-population plus forte pour IT_06 (σ = 0,0186) que pour IT_07 (σ = 0,0089).

**Biologie** : Le delta plus élevé chez IT_07 (altitude 1905 m) suggère que cette population
possède davantage de sites CpG avec des pourcentages intermédiaires proches du seuil 50 %.
En haute altitude, les conditions fluctuantes pourraient induire une méthylation plastique
générant des valeurs modérées plutôt qu'extrêmes. L'approche quantitative est particulièrement
recommandée pour étudier des populations en environnement contraignant ou hétérogène.

---

### Delta par individu

![delta individu](results/figures/delta_par_individu.png)

---

## Validation technique

- **Reproductibilité** : tous les scripts sont documentés et paramétrables. Les fichiers
  intermédiaires sont en `.parquet` (lecture/écriture rapide, conservation du typage).
- **Parallélisation** : l'étape d'intersection utilise `ProcessPoolExecutor` pour traiter
  les 28 paires de fichiers en parallèle.
- **Tests de sensibilité** : les seuils de couverture 5, 10, 15 et 20 ont été testés ;
  les conclusions restent stables.

---

## Licence et contact

Ce travail est sous licence MIT — voir [LICENSE](../LICENSE) à la racine du dépôt.
Contact : Yaya Touré — yaya.toure@etudiant.univ-perp.fr

*Dernière mise à jour : mai 2026*
