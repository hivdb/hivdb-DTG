all:
	./env/bin/python3 work.py

init:
	python -m venv env
	./env/bin/pip3 install -r requirements.txt

ipython:
	./env/bin/ipython

.PHONY: all python init
