### `include/`

This is the root of the include folder, which is where most of the library is 
contained. In addition to this README, this folder holds:

-  `El/`: the supporting header files
-  `El.hpp`: the catch-all header file (users should start by just 
   including this)
-  `El-lite.hpp`: the minimal header file, which can be used in 
   combination with manually including header files from `El/` in order
   to decrease build times
