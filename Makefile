# @file
# @copyright 2025 Terry Golubiewski, all rights reserved.
# @author Terry Golubiewski

export PROJNAME:=ekman
export PROJDIR := $(abspath .)
# SWDEV := $(PROJDIR)/SwDev
include $(SWDEV)/project.mk

APP  :=$(abspath $(HOME)/App)
BOOST:=$(APP)/boost
UNITS:=$(APP)/mp-units
GSL  :=$(APP)/gsl-lite

# Must use "=" instead of ":=" because $E will be defined below.
TJG_E=tjg.$E
LEI_E=lei.$E
TJG2_E=tjg2.$E
TJG3_E=tjg3.$E

TGT1=$(TJG_E)
TGT2=$(LEI_E)
TGT3=$(TJG2_E)
TGT4=$(TJG3_E)
TARGETS=$(TGT1) $(TGT2) $(TGT3) $(TGT4)

SRC1:=tjg.cpp Resample.cpp Smooth.cpp
SRC2:=lei.cpp
SRC3:=tjg2.cpp BoundarySwaths.cpp
SRC4:=tjg3.cpp BoundarySwaths.cpp

SOURCE:=$(SRC1) $(SRC2) $(SRC3) $(SRC4)

SYSINCL:=$(BOOST) $(addsuffix /include, $(UNITS)/src/core $(GSL))
INCLUDE:=$(PROJDIR)

# Must use "=" because LIBS will be changed below.
LIBS=

#DEBUG=1

SHELL:=/bin/bash
.SHELLFLAGS:=-eu -o pipefail -c
.ONESHELL:

include $(SWDEV)/$(COMPILER).mk
include $(SWDEV)/build.mk

.PHONY: all clean scour

all: depend $(TARGETS)

$(TGT1): $(OBJ1) $(LIBS)
        $(LINK)

$(TGT2): $(OBJ2) $(LIBS)
        $(LINK)

$(TGT3): $(OBJ3) $(LIBS)
        $(LINK)

$(TGT4): $(OBJ4) $(LIBS)
        $(LINK)
