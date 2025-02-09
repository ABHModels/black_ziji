# Compiler
CPP := g++

# Flags
GSL_FLAGS := -lgsl -lgslcblas -lm
LD_FLAGS := -fopenmp -L/usr/local/lib/ # Adjust based on system
INC_FLAGS := -I/usr/local/include/ # Adjust based on system
INC_FLAGS += -Izijiray/include

# Add xtensor include path
XTENSOR_SRC_PATH := ./external/xtensor/include
XTL_SRC_PATH := ./external/xtl/include
INC_FLAGS += -I$(XTENSOR_SRC_PATH) -I$(XTL_SRC_PATH)

CPPFLAGS := $(INC_FLAGS) $(LD_FLAGS) $(GSL_FLAGS) -O3 -Wall

# Source files
SRC_RAY := zijiray/src/utils/*.cc zijiray/src/GR/*.cc zijiray/src/raytrace/*.cc

TRANSIT_SRC := $(SRC_RAY) zijiray/src/transit/main.cc
TRANSIT_BIN := $(CPP) $(CPPFLAGS) -shared -fPIC -o bin/ray/libtransit_ray.so $(TRANSIT_SRC)

LPG_SRC := $(SRC_RAY) zijiray/src/lpg/main.cc
LPG_BIN := $(CPP) $(CPPFLAGS) -shared -fPIC -o bin/ray/liblpg.so $(LPG_SRC)

# Targets
.PHONY: all clean lpg_bin transit_bin line_bin reb_bin

lpg_bin:
	$(LPG_BIN)

transit_bin:
	$(TRANSIT_BIN)

SET_LINE_PY := conv_core/core/setup_ray.py
SET_REB_PY := conv_core/core/setup_rebin.py

line_bin:
	python3 $(SET_LINE_PY) build_ext --build-lib bin/conv/
	rm -rf build

reb_bin:
	python3 $(SET_REB_PY) build_ext --build-lib bin/conv/
	rm -rf build

all: lpg_bin transit_bin line_bin reb_bin

clean:
	rm -rf bin/ray/*
	rm -rf bin/conv/*