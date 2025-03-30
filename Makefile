#!/usr/bin/make -f
# Makefile for DISTRHO Plugins #
# ---------------------------- #
# Created by falkTX
#

# --------------------------------------------------------------
# Project name, used for binaries

NAME = cairo_ui
CXXFLAGS += -I./dgl/src/pugl-upstream/include

# --------------------------------------------------------------
# Files to build

FILES_DSP = \
	Fildes.cpp

FILES_UI  = \
	Artwork.cpp \
	FildesUI.cpp

# --------------------------------------------------------------
# Do some magic

UI_TYPE = cairo
include ../../Makefile.plugins.mk
CXXFLAGS += -I/usr/include/pugl-0

# --------------------------------------------------------------
# Enable all possible plugin types

TARGETS = vst2

all: $(TARGETS)

# --------------------------------------------------------------
