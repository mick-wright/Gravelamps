# Construct the lensing source codes into programs to generate lensing data for the interpolant
#
# Targets:
# 	all - Makes lens programs in both current and binary directories
# 	bin_only (default) - Makes lens programs in an editable directory on path
# 	current_only - Makes lens programs in current directory only
# 
# Makefile Written by Mick Wright 2021

# Looking for a compiler being set in the GXX environment variable. If this is not set, use g++
ifeq ($(strip $(GXX)),)
	GXX = g++
endif

# Setting the flags properly. On a local install to default python, flags as set should function.
# In a conda environment additional steps need to be taken to make sure the locations are properly
# done
ifeq ($(strip $(CONDA_DEFAULT_ENV)),)
	GXX_FLAGS = -larb -lflint -O3 -fopenmp
else
	GXX_FLAGS = -isystem $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib -larb -lflint -O3 -fopenmp
endif

# File will attempt to find all .cc files in the folder and create corresponding programs with the same name
CC_FILES = $(wildcard *.cc)
OUT_FILES = $(CC_FILES:%.cc=%) 

CXX_FILES = $(wildcard *.cxx) 
OUT_LIB = $(CXX_FILES:%.cxx=%.so)
# For the binary directory, the makefile will default to using $HOME/bin. It will alert the user if this is not on PATH
BIN_DIR = $(HOME)/bin
LIB_DIR = $(HOME)/.local/lib
OUT_BIN = $(OUT_FILES:%=$(BIN_DIR)/%)
OUT_LIB_FULL = $(OUT_LIB:%=$(LIB_DIR)/%) 
PATH_CHECK = $(findstring $(BIN_DIR),$(PATH))

# Make lens programs in binary directory
bin_only: $(OUT_LIB_FULL) $(OUT_BIN) 

ifeq ($(BIN_DIR),$(PATH_CHECK))
	@echo $(BIN_DIR) on PATH
else
	@echo $(BIN_DIR) not on PATH, will cause issues!
endif

# Make lens programs in current directory only
current_only: $(OUT_LIB_FULL) $(OUT_FILES)

# Make lens programs in both current and binary directories
all: bin_only current_only

# Rule for compiling the .cc source files
%: %.cc
	$(GXX) -o $@ $< $(GXX_FLAGS) -L$(LIB_DIR) -l$(<:%lens.cc=%) -Wl,-rpath=$(LIB_DIR) 

$(BIN_DIR)/%: %.cc
	$(GXX) -o $@ $< $(GXX_FLAGS) -L$(LIB_DIR) -l$(<:%lens.cc=%) -Wl,-rpath=$(LIB_DIR)

$(LIB_DIR)/%.so:%.cxx
	$(GXX) -o $@ $< $(GXX_FLAGS) -fPIC -shared

# Debugging, allows the printing of the variables
variables:
	@echo CC_FILES: $(CC_FILES)
	@echo OUT_FILES: $(OUT_FILES)
	@echo BIN_DIR: $(BIN_DIR)
	@echo OUT_BIN: $(OUT_BIN)
	@echo OUT_LIB_FULL: $(OUT_LIB_FULL)
	@echo LIB_DIR: $(LIB_DIR) 

#CLEAN TARGET
.PHONY: clean
clean:
	rm -f $(OUT_FILES) $(OUT_BIN) $(OUT_LIB_FULL) 

.PHONY: clean_local
clean_local:
	rm -f $(OUT_FILES) $(OUT_LIB_FULL)
