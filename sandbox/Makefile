include ../conf/ElVars

OS := $(shell uname)

all: test
ifeq ($(OS), Linux)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(EL_LIB) ./test
endif
ifeq ($(OS), Darwin )
	DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH):$(EL_LIB) ./test
endif

%: %.cpp
	$(CXX) $(EL_COMPILE_FLAGS) $< -o $@ $(EL_LINK_FLAGS) $(EL_LIBS)

clean:
	rm test

