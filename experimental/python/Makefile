include /home/poulson/Install/conf/elemvariables

libElem.so: interface.cpp
	${CXX} ${ELEM_COMPILE_FLAGS} -c -fPIC interface.cpp -o interface.o
	${CXX} -shared -Wl,-soname,libElem.so -o libElem.so interface.o ${ELEM_LINK_FLAGS} ${ELEM_LIBS}

clean:
	rm -f libElem.so interface.o
