################################
# project configuration
################################

libname = spline

# modules -- header-only
modules_h =

# modules -- header-plus-object 
modules_ho = spline wavefunction_basis wavefunction_class

# programs
##programs = 
programs = h2utils_input

CC := $(CXX)

################################
# common definitions
################################

COMMON_MAKE_DIR ?= .
include $(COMMON_MAKE_DIR)/common.mk

################################
# special dependencies
################################

# program linking
CC := $(CXX)

# external libraries
LDLIBS += -lgsl

CPPFLAGS += -DSPLINE_NO_FANCY_INTEGRATION
CXXFLAGS += -std=c++11
