all:
	./env/bin/python3 work.py

init:
	python -m venv env

ipython:
	./env/bin/ipython

.PHONY: all python init
