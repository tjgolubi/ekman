# @file
# @copyright 2025 Terry Golubiewski, all rights reserved.
# @author Terry Golubiewski

export PROJNAME:=ekman
export PROJDIR := $(abspath ..)
# SWDEV := $(PROJDIR)/SwDev
include $(SWDEV)/project.mk

APP  :=$(abspath $(HOME)/App)
BOOST:=$(APP)/boost
UNITS:=$(APP)/mp-units
GSL  :=$(APP)/gsl-lite
XML  :=$(APP)/pugixml

# Must use "=" instead of ":=" because $E will be defined below.
FIELD_XML_E=field_xml.$E
SWATH_XML_E=swath_xml.$E
TJG_XML_E=tjg_xml.$E
TJG_E=tjg.$E

TGT1=$(FIELD_XML_E)
TGT2=$(SWATH_XML_E)
TGT3=$(TJG_XML_E)
TGT4=$(TJG_E)
TARGETS=$(TGT1) $(TGT2) $(TGT3) $(TGT4)

SRC1:=field_xml.cpp
SRC2:=swath_xml.cpp
SRC3:=tjg_xml.cpp
SRC4:=tjg.cpp
SOURCE:=$(SRC1) $(SRC2) $(SRC3) $(SRC4)

SYSINCL:=$(XML)/src $(BOOST) $(addsuffix /include, $(UNITS)/core $(UNITS)/systems $(GSL))
INCLUDE:=$(PROJDIR)

# Must use "=" because $A Will be defined below.
LIBS=$(XML)/build/libpugixml.$A

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
