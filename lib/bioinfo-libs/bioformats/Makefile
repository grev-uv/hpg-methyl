
# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../..
endif

BIOFORMATS_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats

.PHONY:	all fastq bam-sam vcf gff ped region

all: fastq bam-sam vcf gff ped region

fastq:
	(cd $(BIOFORMATS_DIR)/fastq && make LIBS_HOME=$(LIBS_HOME))

bam-sam:
	cd $(BIOFORMATS_DIR)/bam-sam && make LIBS_HOME=$(LIBS_HOME)

vcf:
	cd $(BIOFORMATS_DIR)/vcf && make LIBS_HOME=$(LIBS_HOME)

ped:
	cd $(BIOFORMATS_DIR)/ped && make LIBS_HOME=$(LIBS_HOME)

gff:
	cd $(BIOFORMATS_DIR)/gff && make LIBS_HOME=$(LIBS_HOME)

family:
	cd $(BIOFORMATS_DIR)/family && make LIBS_HOME=$(LIBS_HOME)

region:
	cd $(BIOFORMATS_DIR)/features/region && make LIBS_HOME=$(LIBS_HOME)


