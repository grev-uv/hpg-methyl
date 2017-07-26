
# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../..
endif

ALIGNERS_DIR = $(LIBS_ROOT)/bioinfo-libs/aligners

.PHONY:	all bwt sw

all: bwt sw

bwt:
	cd $(ALIGNERS_DIR)/bwt && make LIBS_HOME=$(LIBS_HOME)

sw:
	cd $(ALIGNERS_DIR)/sw && make LIBS_HOME=$(LIBS_HOME)

