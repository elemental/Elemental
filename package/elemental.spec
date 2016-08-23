Name:	elemental
Version:	0.86
Release:	1%{?dist}
Summary:	Library for distributed-memory dense/sparse-direct linear algebra 
Group:	Development/Libraries
License:	BSD
URL:	http://libelemental.org
Source0:	https://github.com/rhl-/Elemental/archive/%{version}-rc4.tar.gz 

BuildRequires: cmake
BuildRequires: metis-devel >= 5.1.0
BuildRequires: openblas-devel
BuildRequires: python2-devel 
BuildRequires: qd-devel

%{?el6:BuildRequires:  devtoolset-4-toolchain}
%{?el7:BuildRequires:  devtoolset-4-toolchain}

%description 
A modern C++ library for distributed-memory linear algebra.

%package common
Summary: Files in common between mpich and openmpi
Group: Development/Libraries
%description common 
Files not specific to mpich or openmpi

%package devel 
Summary: Elemental C/C++ Header Files
Group: Development/Libraries
%description devel
Use this package for building off of Elemental

%package python2 
Summary: Python 2 Bindings 
Group: Development/Libraries
%description python2
This package contains the python bindings for using Elemental through a python shell

%package openmpi
Summary: OpenMPI variant of Elemental
Group: Development/Libraries
BuildRequires: openmpi-devel
BuildRequires: scalapack-openmpi-devel
# Require explicitly for dir ownership and to guarantee the pickup of the right runtime
Requires: openmpi
Requires: %{name}-common = %{version}-%{release}
%description openmpi
Contains the library, unit tests, and example drivers built against OpenMPI

%package mpich
Summary: MPICH variant of Elemental
Group: Development/Libraries
BuildRequires: mpich-devel
BuildRequires: scalapack-mpich-devel
# Require explicitly for dir ownership and to guarantee the pickup of the right runtime
Requires: mpich
Requires: %{name}-common = %{version}-%{release}
%description mpich
Contains the library, unit tests, and example drivers built against MPICH

%prep
%autosetup 

%build

%if 0%{?rhel}
source /opt/rh/devtoolset-4/enable
%endif

%define dobuild() \
mkdir $MPI_COMPILER; \
cd $MPI_COMPILER;  \
%cmake -DCMAKE_BUILD_TYPE=Release -DNO_BINARY_SUBDIRECTORIES=True -DCMAKE_RELEASE_POSTFIX="$MPI_SUFFIX" -DCMAKE_EXECUTABLE_SUFFIX_CXX="$MPI_SUFFIX" -DEL_TESTS=ON -DEL_EXAMPLES=ON -DINSTALL_PYTHON_PACKAGE=ON -DGFORTRAN_LIB="$(gfortran -print-file-name=libgfortran.so)" -DEL_DISABLE_SCALAPACK=ON -DEL_DISABLE_PARMETIS=ON .. ; \
make %{?_smp_mflags}; \
cd .. ; \


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
%{_prefix}/%{_sysconfdir}/elemental/CMake/*

%files python2
%{python2_sitelib}/*

# All files shared between the serial and different MPI versions
%files common 
%{_datadir}/elemental/*
%{_datadir}/doc/Elemental/*

# All openmpi linked files
%files openmpi 
%{_libdir}/*_openmpi.*
%{_bindir}/*_openmpi*

# All mpich linked files
%files mpich 
%{_libdir}/*_mpich.*
%{_bindir}/*_mpich*

%changelog
* Thu Jul 28 2016 Ryan H. Lewis <me@ryanlewis.net> - 0.86-1
- Initial RPM
