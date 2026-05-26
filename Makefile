PYTHON  = venv/bin/python3
SCRIPTS = scripts
RESULTS = results

.PHONY: all intersection filtrage h2_qual h2_quant figures clean help install

help:
	@echo "============================================"
	@echo "  Pipeline Heritabilite Methylation CpG"
	@echo "============================================"
	@echo "  make install       -> Installer les deps"
	@echo "  make intersection  -> Etape 1 : P0 inter R2"
	@echo "  make filtrage      -> Etape 2 : Filtrage coverage"
	@echo "  make h2_qual       -> Etape 3A : h² qualitatif"
	@echo "  make h2_quant      -> Etape 3B : h² quantitatif"
	@echo "  make figures       -> Etape 4 : Visualisation"
	@echo "  make all           -> Toutes les etapes"
	@echo "  make clean         -> Supprimer resultats"
	@echo "============================================"

install:
	pip install polars pyarrow tqdm pandas matplotlib seaborn scipy

all: intersection filtrage h2_qual h2_quant figures

intersection:
	@echo ">>> Etape 1 : Intersection P0/R2..."
	$(PYTHON) $(SCRIPTS)/01_intersection_CpG.py

filtrage:
	@echo ">>> Etape 2 : Filtrage couverture..."
	$(PYTHON) $(SCRIPTS)/02_filtrage_couverture.py

h2_qual:
	@echo ">>> Etape 3A : Heritabilite qualitative..."
	$(PYTHON) $(SCRIPTS)/03A_heritabilite_qualitative.py

h2_quant:
	@echo ">>> Etape 3B : Heritabilite quantitative..."
	$(PYTHON) $(SCRIPTS)/03B_heritabilite_quantitative.py

figures:
	@echo ">>> Etape 4 : Visualisation..."
	$(PYTHON) $(SCRIPTS)/04_visualisation.py

clean:
	rm -f $(RESULTS)/*.parquet $(RESULTS)/*.csv
	rm -f $(RESULTS)/filtrage/*.parquet $(RESULTS)/filtrage/*.csv
	rm -f $(RESULTS)/heritabilite/*.csv
	rm -f $(RESULTS)/figures/*.png
	@echo "Resultats supprimes."
