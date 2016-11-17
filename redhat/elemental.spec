Name:	elemental
Version:	0.87
Release:	2%{?dist}
Summary:	Library for distributed-memory dense/sparse-direct linear algebra 
Group:	Development/Libraries
License:	BSD
URL:	http://libelemental.org
Source0:	https://github.com/elemental/Elemental/archive/master.zip 

BuildRequires: cmake
BuildRequires: metis-devel >= 5.1.0
BuildRequires: openblas-devel
BuildRequires: python2-devel 
BuildRequires: qd-devel
BuildRequires: qt5-qtbase-devel
BuildRequires: gmp-devel

%{?el6:BuildRequires:  devtoolset-4-toolchain}
%{?el7:BuildRequires:  devtoolset-4-toolchain}

%description 
A modern C++ library for distributed-memory linear algebra.

%package common
Summary: Files in common between mpich and openmpi
Group: Development/Libraries
Requires: qt5-qtbase
%description common 
Files not specific to mpich or openmpi

%package devel 
Summary: Elemental C/C++ Header Files
Group: Development/Libraries
%description devel
Use this package for building off of Elemental

## OpenMPI Subpackages
%package openmpi
Summary: OpenMPI variant of Elemental
Group: Development/Libraries
BuildRequires: openmpi-devel
Requires: openmpi
Requires: %{name}-common = %{version}-%{release}
%description openmpi
Contains the library, built against OpenMPI

%package openmpi-examples
Summary: OpenMPI variant of Elemental
Group: Development/Libraries
Requires: %{name}-openmpi = %{version}-%{release}
%description openmpi-examples
Contains the example drivers built against OpenMPI

%package openmpi-tests
Summary: OpenMPI variant of Elemental
Group: Development/Libraries
Requires: %{name}-openmpi = %{version}-%{release}
%description openmpi-tests
Contains the test drivers built against OpenMPI

%package -n python2-elemental-openmpi 
Summary: Python 2 Bindings 
Group: Development/Libraries
Requires: %{name}-openmpi = %{version}-%{release}
%description -n python2-elemental-openmpi
This package contains the python bindings for using Elemental through a python shell with OpenMPI

## MPICH Subpackages
%package mpich
Summary: MPICH variant of Elemental
Group: Development/Libraries
BuildRequires: mpich-devel
Requires: mpich
Requires: %{name}-common = %{version}-%{release}
%description mpich
Contains the library, unit tests, and example drivers built against MPICH

%package mpich-examples
Summary: MPICH variant of Elemental
Group: Development/Libraries
Requires: %{name}-mpich = %{version}-%{release}
%description mpich-examples
Contains the example drivers built against MPICH

%package mpich-tests
Summary: MPICH variant of Elemental
Group: Development/Libraries
Requires: %{name}-mpich = %{version}-%{release}
%description mpich-tests
Contains the test drivers built against MPICH

%package -n python2-elemental-mpich
Summary: Python 2 Bindings 
Group: Development/Libraries
Requires: %{name}-mpich = %{version}-%{release}
%description -n python2-elemental-mpich
This package contains the python bindings for using Elemental through a python shell with MPICH

%prep
%autosetup 

%build

%if 0%{?rhel}
source /opt/rh/devtoolset-4/enable
%endif

%define dobuild() \
mkdir $MPI_COMPILER; \
cd $MPI_COMPILER;  \
%cmake -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpic++" -DCMAKE_BUILD_TYPE=Release -DBUILD_METIS=OFF -DEL_USE_QT5=ON -DBINARY_SUBDIRECTORIES=False -DEL_TESTS=ON -DEL_EXAMPLES=ON -DINSTALL_PYTHON_PACKAGE=ON -DGFORTRAN_LIB="$(gfortran -print-file-name=libgfortran.so)" -DEL_DISABLE_PARMETIS=ON -DCMAKE_INSTALL_BINDIR="$MPI_BIN" -DCMAKE_INSTALL_LIBDIR="$MPI_LIB" -DPYTHON_SITE_PACKAGES="$MPI_PYTHON_SITEARCH" .. ; \
make %{?_smp_mflags}; \
ls
find . -name "libEl.so.0"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.; \
env CTEST_OUTPUT_ON_FAILURE=1 make test \
cd .. ; \

# Set compiler variables to MPI wrappers
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

#The Elemental headers
%files devel
%{_includedir}/*
%{_prefix}/%{_sysconfdir}/elemental/CMake/*

# All files shared between the serial and different MPI versions
%files common 
%{_datadir}/elemental/*
%{_datadir}/doc/Elemental/*

# All openmpi linked files
%files openmpi 
%{_libdir}/openmpi/lib/*

%files openmpi-examples
%{_libdir}/openmpi/bin/examples/*

%files openmpi-tests
%{_libdir}/openmpi/bin/tests/*

%files -n python2-elemental-openmpi
%{python2_sitearch}/openmpi/*

# All mpich files
%files mpich 
%{_libdir}/mpich/lib/*

%files mpich-examples
%{_libdir}/mpich/bin/examples/*

%files mpich-tests
%{_libdir}/mpich/bin/tests/*

%files -n python2-elemental-mpich 
%{python2_sitearch}/mpich/*

%changelog
* Sat Oct 29 2016 Ryan H. Lewis <me@ryanlewis.net> - 0.87-1
- Dropped Scalapack 
- Enabled Qt5
- updated Source0 to master

* Thu Jul 28 2016 Ryan H. Lewis <me@ryanlewis.net> - 0.86-1
- Initial RPM
