### `include/`

This is the root of the include folder, which is where most of the library is 
contained. In addition to this README, this folder holds:

-  `elemental/`: the supporting header files
-  `elemental.hpp`: the catch-all header file (users should start by just 
   including this)
-  `elemental-lite.hpp`: the minimal header file, which can be used in 
   combination with manually including header files from `elemental/` in order
   to decrease build times
