Name:           elemental
Version:	0.86
Release:        1%{?dist}
Summary: Elemental is an open-source library for distributed-memory dense and sparse-direct linear algebra 
Group: Development/Libraries
License: BSD        
URL: http://libelemental.org           
Source0: https://github.com/rhl-/Elemental/archive/%{version}-rc4.tar.gz 

%if 0%{?rhel} ==5 
%global package_cmake %cmake28
BuildRequires:  cmake28
%endif
%if 0%{?rhel} !=5 
%global package_cmake %cmake
BuildRequires: cmake
%endif

BuildRequires: metis-devel >= 5.1.0
BuildRequires: openblas-devel
BuildRequires: python2-devel 
BuildRequires: qd-devel

%{?el5:BuildRequires:  devtoolset-2}
%{?el6:BuildRequires:  devtoolset-4}
%{?el7:BuildRequires:  devtoolset-4}

%package common
Summary: common stuff between both elemental versions
Group: Development/Libraries
%description common

%package devel 
Summary: devel headers
Group: Development/Libraries
%description devel

%package python2 
Summary: python bindings 
Group: Development/Libraries
%description python2

%package openmpi
Summary: openmpi elemental 
Group: Development/Libraries
BuildRequires: openmpi-devel
BuildRequires: scalapack-openmpi-devel
# Require explicitly for dir ownership and to guarantee the pickup of the right runtime
Requires: openmpi
Requires: %{name}-common = %{version}-%{release}
%description openmpi

%package mpich
Summary: mpich elemental 
Group: Development/Libraries
BuildRequires: mpich-devel
BuildRequires: scalapack-mpich-devel
# Require explicitly for dir ownership and to guarantee the pickup of the right runtime
Requires: mpich
Requires: %{name}-common = %{version}-%{release}
%description mpich


%description
Elemental is an open-source library for distributed-memory dense and sparse-direct linear algebra and optimization which builds on top of BLAS, LAPACK, and MPI using modern C++ and additionally exposes interfaces to C and Python (with a Julia interface beginning development). The development of Elemental has led to a number of research articles and a number of related projects, such as the parallel sweeping preconditioner, PSP, and a parallel algorithm for Low-rank Plus Sparse MRI, RT-LPS-MRI.

%prep
%autosetup 

%build

%if 0%{?rhel} > 5
scl enable devtoolset-4 bash
%endif
%if 0%{?rhel} == 5
scl enable devtoolset-2 bash
%endif

%define dobuild() \
mkdir $MPI_COMPILER; \
cd $MPI_COMPILER;  \
%package_cmake -DCMAKE_BUILD_TYPE=Release -DNO_BINARY_SUBDIRECTORIES=True -DCMAKE_RELEASE_POSTFIX="$MPI_SUFFIX" -DCMAKE_EXECUTABLE_SUFFIX_CXX="$MPI_SUFFIX" -DEL_TESTS=ON -DEL_EXAMPLES=ON -DINSTALL_PYTHON_PACKAGE=ON -DGFORTRAN_LIB="$(gfortran -print-file-name=libgfortran.so)" -DEL_DISABLE_SCALAPACK=ON -DEL_DISABLE_PARMETIS=ON .. ; \
make %{?_smp_mflags}; \
cd .. ; \

# Build parallel versions: set compiler variables to MPI wrappers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77

## Build OpenMPI version
%{_openmpi_load}
%dobuild
%{_openmpi_unload}

# Build mpich version
%{_mpich_load}
%dobuild
%{_mpich_unload}


%install
## Install OpenMPI version
%{_openmpi_load}
make -C $MPI_COMPILER install/fast DESTDIR=%{buildroot} INSTALL="install -p" CPPROG="cp -p"
%{_openmpi_unload}

# Install MPICH2 version
%{_mpich_load}
make -C $MPI_COMPILER install/fast DESTDIR=%{buildroot} INSTALL="install -p" CPPROG="cp -p"
%{_mpich_unload}

rm -rf %{buildroot}/%{_prefix}/conf

%files devel
%{_includedir}/*

%files python2
%{python2_sitelib}/*

# All files shared between the serial and different MPI versions
%files common 
%{_datadir}/*
%{_prefix}/%{_sysconfdir}/elemental/CMake/*

# All openmpi linked files
%files openmpi 
%{_libdir}/*_openmpi.*
%{_bindir}/*_openmpi*

# All mpich linked files
%files mpich 
%{_libdir}/*_mpich.*
%{_bindir}/*_mpich*

%changelog
* Thu Jul 28 2016 Initial RPM
- 
