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

TGT1=$(TJG_E)
TGT2=$(LEI_E)
TARGETS=$(TGT1) $(TGT2)

SRC1:=tjg.cpp Resample.cpp Smooth.cpp
SRC2:=lei.cpp

SOURCE:=$(SRC1) $(SRC2)

SYSINCL:=$(BOOST) $(addsuffix /include, $(UNITS) $(GSL))
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
