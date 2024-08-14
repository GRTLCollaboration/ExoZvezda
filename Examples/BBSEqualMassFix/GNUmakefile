# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

#ifndef BOOST_ROOT
#		$(error The BosonStar example requires BOOST_ROOT to be set to the root of the boost library directory.)
#endif

GRCHOMBO_SOURCE = ../../Source

ebase := Main_BosonStar

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/Matter  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/BosonStars  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BosonStars \
            $(GRCHOMBO_SOURCE)/InitialConditions/BosonStars/EqualMassFix \
            $(GRCHOMBO_SOURCE)/ApparentHorizonFinder 

include $(CHOMBO_HOME)/mk/Make.test
