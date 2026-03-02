# Makefile for hungarian-algorithm-fortran
FC = gfortran
FFLAGS = -O3 -Wall -Wextra

# Platform-specific shared library extension and flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIB_EXT = dylib
	SHARED_FLAG = -shared -fPIC -march=native
else ifeq ($(UNAME_S),Linux)
	LIB_EXT = so
	SHARED_FLAG = -shared -fPIC -march=native
else
	# Windows / MinGW
	LIB_EXT = dll
	SHARED_FLAG = -shared
endif

LIB_NAME = libhungarian.$(LIB_EXT)
TEST_BIN = test_exe

.PHONY: all lib test clean

all: lib

lib: $(LIB_NAME)

$(LIB_NAME): hungarian_algorithm.f90
	$(FC) $(SHARED_FLAG) $(FFLAGS) $< -o $@

# Phony target to compile and run the test
test: $(TEST_BIN)
	./$(TEST_BIN)

# Rule to build the actual test executable
$(TEST_BIN): test.f90 hungarian_algorithm.f90
	$(FC) $(FFLAGS) hungarian_algorithm.f90 test.f90 -o $@

clean:
	rm -f *.o *.mod *.smod $(LIB_NAME) $(TEST_BIN) *.exe