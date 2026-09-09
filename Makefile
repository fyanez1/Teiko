# Loblaw Bio immune cell-count analysis
#
#   make setup      install Python dependencies
#   make pipeline   build cell_counts.db from cell-count.csv and run Parts 2-4 (writes outputs/)
#   make dashboard  start the interactive Streamlit dashboard (http://localhost:8501)
#   make test       run the automated checks
#   make clean      remove generated files

PYTHON ?= python3
STREAMLIT_PORT ?= 8501

.PHONY: setup pipeline dashboard test clean

setup:
	PIP_BREAK_SYSTEM_PACKAGES=1 $(PYTHON) -m pip install -r requirements.txt || \
	PIP_BREAK_SYSTEM_PACKAGES=1 $(PYTHON) -m pip install --user -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard/app.py --server.port $(STREAMLIT_PORT) --server.headless true

test:
	$(PYTHON) -m pytest -q

clean:
	rm -f cell_counts.db
	rm -rf outputs/* .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
