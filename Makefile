all:
	./env/bin/python3 work.py

init:
	python3 -m venv env
	@echo "env" >> .gitignore
	./env/bin/pip3 install -r requirements.txt

freeze:
	@./env/bin/pip3 freeze > requirements.txt

ipython:
	./env/bin/ipython

.PHONY: all init freeze ipython
