# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = .

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# api doc
api:
	@echo "Generating API documentation..."
	@sphinx-apidoc -o source ../pyCADD

# Distributed build targets
distributed:
	@echo "Building documentation with distributed environments..."
	@python3 build_docs.py

clean-distributed:
	@echo "Cleaning distributed build artifacts..."
	@rm -rf build/

check-env:
	@echo "Checking environments..."
	@python3 build_docs.py --check

# Build individual modules
build-dock:
	@echo "Building Dock module documentation..."
	@python3 build_docs.py --module Dock

build-dynamic:
	@echo "Building Dynamic module documentation..."
	@python3 build_docs.py --module Dynamic

build-dance:
	@echo "Building Dance module documentation..."  
	@python3 build_docs.py --module Dance

build-demand:
	@echo "Building Demand module documentation..."
	@python3 build_docs.py --module Demand

build-density:
	@echo "Building Density module documentation..."
	@python3 build_docs.py --module Density

.PHONY: help Makefile distributed clean-distributed build-dock build-dynamic

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
