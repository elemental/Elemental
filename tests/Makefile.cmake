include @CMAKE_INSTALL_PREFIX@/conf/ElVars

OS := $(shell uname)

TESTS := $(shell ls *.cpp | sed s/\.cpp//g)

.PHONY: all clean

all: $(TESTS)

%: %.cpp
	$(CXX) $(EL_COMPILE_FLAGS) $< -o $@ $(EL_LINK_FLAGS) $(EL_LIBS)
	if test x$(OS) = xLinux; then LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(EL_LIB) ./$@; fi
	if test x$(OS) = xDarwin; then DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH):$(EL_LIB) ./$@; fi
	rm $@

clean:
	rm -f $(TESTS)

