# Construct the lensing model source codes into a library of functions to be used by the python
# front end of Gravelamps. 
#
# Targets:
# 	all - makes all lensing libraries
#
# Makefile written by Mick Wright 2022

# Compiler checking - if GXX is set use, otherwise set to g++
ifeq ($(strip $(GXX)),)
	GXX = g++
endif

# The default list of flags to set, adjusted depending on whether the makefile is run with or
# outwith a conda environment. Flags set are the inclusion of the flint and arb libraries,
# as many optimisations as possible, and multiprocessing via OpenMP
ifeq ($(strip $CONDA_DEFAULT_ENV)),)
	GXX_FLAGS = -larb -lflint -O3 -fopenmp -g
else
	GXX_FLAGS = -isystem $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib -larb -lflint\
		    -O3 -fopenmp -g
endif

# Libraries to create are determined by the headers in the src directory
HEADER_FILES = $(filter-out src/utils.h, $(wildcard src/*.h))
MODELS = $(HEADER_FILES:src/%.h=lib%)
LIBRARY_NAMES = $(MODELS:%=%.so)
LIBRARY_DIRECTORY = lib
LIBRARY_TARGETS = $(LIBRARY_NAMES:%=$(LIBRARY_DIRECTORY)/%) 

# Rules for creating libraries
$(LIBRARY_DIRECTORY)/%.so :
	@echo '**********************' 
	$(eval LIBRARY_NAME = $(@:$(LIBRARY_DIRECTORY)/%.so=%))
	@echo BUILDING $(LIBRARY_NAME)
	$(eval SOURCE_FILES = $(wildcard $(LIBRARY_NAME)_*.cxx)) 
	@echo SOURCE_FILES : $(SOURCE_FILES)
	$(GXX) -o $@ $(SOURCE_FILES) utils.cxx $(GXX_FLAGS) -fPIC -shared
	@echo $(LIBRARY_NAME) BUILT 
	@echo '**********************'

all: $(LIBRARY_TARGETS)

# Debugging, print all variables
variables:
	@echo HEADER_FILES: $(HEADER_FILES)
	@echo MODELS: $(MODELS)
	@echo LIBRARY_NAMES: $(LIBRARY_NAMES)
	@echo LIBRARY_DIRECTORY: $(LIBRARY_DIRECTORY)
	@echo $(LIBRARY_DIRECTORY)/$(LIBRARY_NAMES)

#CLEAN TARGET
.PHONY: clean
clean:
	rm -f $(LIBRARY_TARGETS)
