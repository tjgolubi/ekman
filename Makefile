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

TGT1=$(TJG_E)
TARGETS=$(TGT1)

SRC1:=tjg.cpp Resample.cpp Smooth.cpp

SOURCE:=$(SRC1)

SYSINCL:=$(BOOST) $(addsuffix /include, $(UNITS)/core $(UNITS)/systems $(GSL))
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
