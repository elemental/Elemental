include /home/poulson/Install/conf/elemvariables

libElem.so: interface.cpp
	${CXX} ${ELEM_COMPILE_FLAGS} -c -fPIC interface.cpp -o elem.o
	${CXX} -shared -Wl,-soname,libElem.so -o libElem.so elem.o ${ELEM_LINK_FLAGS} ${ELEM_LIBS}
