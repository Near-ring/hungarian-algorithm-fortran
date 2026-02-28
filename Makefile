# Makefile for hungarian-algorithm-fortran
FC = gfortran
FFLAGS = -O3 -march=native -fPIC -Wall -Wextra

# Platform-specific shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIB_EXT = dylib
    SHARED_FLAG = -dynamiclib
else ifeq ($(UNAME_S),Linux)
    LIB_EXT = so
    SHARED_FLAG = -shared
else
    LIB_EXT = dll
    SHARED_FLAG = -shared
endif

LIB_NAME = libhungarian.$(LIB_EXT)
TEST_NAME = test_hungarian

.PHONY: all lib test clean

all: lib

lib: $(LIB_NAME)

$(LIB_NAME): hungarian.f90
	$(FC) $(SHARED_FLAG) $(FFLAGS) $< -o $@

# Fortran test program
test: $(TEST_NAME)
	./$(TEST_NAME)

$(TEST_NAME): test_hungarian.f90 hungarian.f90
	$(FC) $(FFLAGS) hungarian.f90 test_hungarian.f90 -o $@

clean:
	rm -f *.o *.mod *.smod $(LIB_NAME) $(TEST_NAME)
