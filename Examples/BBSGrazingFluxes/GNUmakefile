# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

# The user needs to specify GRCHOMBO_SOURCE location, e.g. using 'export GRCHOMBO_SOURCE = /home/folder/GRChombo/Source'

export EXOZVEZDA_SOURCE = $(shell pwd)/../../Source

ebase := Main_BBSGrazingFluxes

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/Matter  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/ApparentHorizonFinder \
            $(EXOZVEZDA_SOURCE)/InitialConditions/BosonStars \
            $(EXOZVEZDA_SOURCE)/InitialConditions/BosonStars/EqualMassFix \
            $(EXOZVEZDA_SOURCE)/Matter \
            $(EXOZVEZDA_SOURCE)/BosonStars \
            $(EXOZVEZDA_SOURCE)/utils \
            $(EXOZVEZDA_SOURCE)/TaggingCriteria 

include $(CHOMBO_HOME)/mk/Make.test
