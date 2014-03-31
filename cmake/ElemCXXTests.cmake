include(FindCXXFeatures)
message(STATUS "CXX11_COMPILER_FLAGS=${CXX11_COMPILER_FLAGS}")
set(EXTRA_FLAGS "${CXX11_COMPILER_FLAGS}")
set(CMAKE_REQUIRED_FLAGS ${CXX11_COMPILER_FLAGS})
set(NORMAL_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::normal_distribution<double> dist(0,1);
         const double x = dist(mt);
         return 0;
     }")
set(UNIFORM_INT_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::uniform_int_distribution<int> dist(0,10);
         const int x = dist(mt);
         return 0;
     }")
set(UNIFORM_REAL_CODE
    "#include <random>
     int main()
     {
         std::random_device rd;
         std::mt19937 mt(rd());
         std::uniform_real_distribution<double> dist(0,1);
         const double x = dist(mt);
         return 0;
     }")
set(STEADYCLOCK_CODE
    "#include <chrono>
     int main()
     {
         std::chrono::steady_clock clock;
         return 0;
     }")
set(ALIAS_CODE
    "#include <complex>
     template<typename Real>
     using Complex = std::complex<Real>;
     int main()
     {
         Complex<double> a;
     }")
set(NOEXCEPT_CODE
    "#include <vector>
     void Foo( const std::vector<int>& x ) noexcept { }
     int main()
     { return 0; }")
check_cxx_source_compiles("${NORMAL_CODE}" ELEM_HAVE_NORMAL_DIST)
check_cxx_source_compiles("${UNIFORM_INT_CODE}" ELEM_HAVE_UNIFORM_INT_DIST)
check_cxx_source_compiles("${UNIFORM_REAL_CODE}" ELEM_HAVE_UNIFORM_REAL_DIST)
check_cxx_source_compiles("${STEADYCLOCK_CODE}" ELEM_HAVE_STEADYCLOCK)
check_cxx_source_compiles("${ALIAS_CODE}" ELEM_HAVE_TEMPLATE_ALIAS)
check_cxx_source_compiles("${NOEXCEPT_CODE}" ELEM_HAVE_NOEXCEPT)
if(NOT ELEM_HAVE_NORMAL_DIST OR NOT ELEM_HAVE_TEMPLATE_ALIAS)
  message(FATAL_ERROR "C++11 random number generation and template aliasing not found. You may want to make sure that your compiler is sufficiently up-to-date (e.g., g++ >= 4.7)")
endif()
# Check for restrict
set(RESTRICT_CODE "int main() { int* RESTRICT a; return 0; }")
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict)
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=restrict")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE_restrict)
if(HAVE___restrict__)
  set(RESTRICT "__restrict__")
  message(STATUS "Using __restrict__ keyword.")
elseif(HAVE___restrict)
  set(RESTRICT "__restrict")
  message(STATUS "Using __restrict keyword.")
elseif(HAVE_restrict)
  set(RESTRICT "restrict")
  message(STATUS "Using restrict keyword.")
else()
  set(RESTRICT "")
  message(STATUS "Could not find a restrict keyword.")
endif()
