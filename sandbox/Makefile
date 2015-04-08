include ../conf/ElVars

OS := $(shell uname)

all: test
	if [ x$(OS) == xLinux ]; then LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(EL_LIB) ./test; fi
	if [ x$(OS) == xDarwin ]; then DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH):$(EL_LIB) ./test; fi

test: test.cpp
	$(CXX) $(EL_COMPILE_FLAGS) $< -o $@ $(EL_LINK_FLAGS) $(EL_LIBS)

clean:
	rm test
