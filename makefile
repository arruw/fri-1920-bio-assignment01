SHELL := /bin/bash
VENV := .env/bin

install:
	python3 -m venv .env
	$(VENV)/pip install -r requirements.txt

run:
	$(VENV)/python3 src/solution.py