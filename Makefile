name := $(shell basename $(PWD))

G4TARGET := $(name)
G4EXLIB := true
G4_NO_VERBOSE := true
CPPFLAGS += $(shell root-config --cflags)
EXTRALIBS := $(shell root-config --ldflags) $(shell root-config --libs)

ifdef VIS
CPPFLAGS += -DVIS
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
