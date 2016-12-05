Name:	elemental
Version:	0.87.5
Release:	2%{?dist}
Summary:	Library for distributed-memory dense/sparse-direct linear algebra 
Group:	Development/Libraries
License:	BSD and Boost and MIT and LGPLv2
URL:	http://libelemental.org
Source0:	https://github.com/elemental/Elemental/archive/v%{version}.tar.gz
#This is excluded to due a compiler bug in PPC:
#gcc.gnu.org/bugzilla/show_bug.cgi?id=78636
ExcludeArch: %{power64}

BuildRequires: environment-modules
BuildRequires: cmake
BuildRequires: metis-devel >= 5.1.0
BuildRequires: openblas-devel
BuildRequires: python2-devel 
BuildRequires: qd-devel
BuildRequires: qt5-qtbase-devel
BuildRequires: gmp-devel
BuildRequires: mpfr-devel
BuildRequires: libmpc-devel

%{?el6:BuildRequires:  devtoolset-4-toolchain}
%{?el7:BuildRequires:  devtoolset-4-toolchain}

%description 
A modern C++ library for distributed-memory linear algebra.

%package common
Summary: Files in common between mpich and openmpi
Group: Development/Libraries
BuildArch: noarch
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

%package openmpi-devel
Summary: OpenMPI variant of Elemental
Group: Development/Libraries
Requires: %{name}-openmpi%{?_isa} = %{version}-%{release}
%description openmpi-devel
Contains the library, built against OpenMPI

%package openmpi-examples
Summary: OpenMPI variant of Elemental
Group: Development/Libraries
Requires: %{name}-openmpi%{?_isa} = %{version}-%{release}
%description openmpi-examples
Contains the example drivers built against OpenMPI

%package -n python2-elemental-openmpi 
Summary: Python 2 Bindings 
Group: Development/Libraries
Requires: %{name}-openmpi%{?_isa} = %{version}-%{release}
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
Contains the library, and example drivers built against MPICH

%package mpich-devel
Summary: MPICH variant of Elemental
Group: Development/Libraries
BuildRequires: mpich-devel
Requires: mpich
Requires: %{name}-mpich%{?_isa} = %{version}-%{release}
%description mpich-devel
Contains the library built against MPICH

%package mpich-examples
Summary: MPICH variant of Elemental
Group: Development/Libraries
Requires: %{name}-mpich%{?_isa} = %{version}-%{release}
%description mpich-examples
Contains the example drivers built against MPICH

%package -n python2-elemental-mpich
Summary: Python 2 Bindings 
Group: Development/Libraries
Requires: %{name}-mpich%{?_isa} = %{version}-%{release}
%description -n python2-elemental-mpich
This package contains the python bindings for using Elemental through a python shell with MPICH

%prep
%autosetup -c -n Elemental-%{version}
mv $(ls -d */|head -n 1)/* . 

%build

%if 0%{?rhel}
source /opt/rh/devtoolset-4/enable
%endif

%define dobuild() \
mkdir $MPI_COMPILER; \
cd $MPI_COMPILER;  \
export CXXFLAGS="%{optflags} -Wl,--as-needed"; \
%cmake -DINSTALL_CMAKE_DIR="%{_libdir}/cmake/" -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpic++" -DCMAKE_BUILD_TYPE=Release -DBUILD_METIS=OFF -DEL_USE_QT5=ON -DBINARY_SUBDIRECTORIES=False -DEL_TESTS=ON -DEL_EXAMPLES=ON -DINSTALL_PYTHON_PACKAGE=ON -DGFORTRAN_LIB="$(gfortran -print-file-name=libgfortran.so)" -DEL_DISABLE_PARMETIS=ON -DCMAKE_INSTALL_BINDIR="$MPI_BIN" -DCMAKE_INSTALL_LIBDIR="$MPI_LIB" -DPYTHON_SITE_PACKAGES="$MPI_PYTHON_SITEARCH" .. ; \
make %{?_smp_mflags}; \
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

%check 

%define docheck() \
export  CTEST_OUTPUT_ON_FAILURE=1; \
cd $MPI_COMPILER ; \
export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH; \
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd):$(pwd)/external/pmrrr:$(pwd)/external/suite_sparse; \
ctest -V %{?_smp_mflags}; \
export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH; \
cd .. ; \

## Build OpenMPI version
%{_openmpi_load}
%docheck
%{_openmpi_unload}

# Build mpich version
%{_mpich_load}
%docheck
%{_mpich_unload}

%install
## Install OpenMPI version
%{_openmpi_load}
make -C $MPI_COMPILER install/fast DESTDIR=%{buildroot} INSTALL="install -p" CPPROG="cp -p"
rm -f %{buildroot}/$MPI_BIN/tests-*
%{_openmpi_unload}

# Install MPICH2 version
%{_mpich_load}
make -C $MPI_COMPILER install/fast DESTDIR=%{buildroot} INSTALL="install -p" CPPROG="cp -p"
rm -f ${buildroot}/$MPI_BIN/tests-*
%{_mpich_unload}

mv %{buildroot}%{_docdir}/Elemental %{buildroot}%_pkgdocdir
rm -rf %{buildroot}/%{_prefix}/conf

#The Elemental headers
%files devel
%{_includedir}/*
%{_libdir}/cmake/elemental/*

# All files shared between the serial and different MPI versions
%files common 
%{_datadir}/elemental/*
%_pkgdocdir/
%license debian/copyright
%license LICENSE

# All openmpi linked files
%files openmpi 
%{_libdir}/openmpi/lib/*.so.*

%files openmpi-devel
%{_libdir}/openmpi/lib/*.so

%files openmpi-examples
%{_libdir}/openmpi/bin/*

%files -n python2-elemental-openmpi
%{python2_sitearch}/openmpi/*

# All mpich files
%files mpich 
%{_libdir}/mpich/lib/*.so.*

%files mpich-devel
%{_libdir}/mpich/lib/*.so

%files mpich-examples
%{_libdir}/mpich/bin/*

%files -n python2-elemental-mpich 
%{python2_sitearch}/mpich/*

%license debian/copyright

%changelog
* Sat Oct 29 2016 Ryan H. Lewis <me@ryanlewis.net> - 0.87-1
- Dropped Scalapack 
- Enabled Qt5
- updated Source0 to master

* Thu Jul 28 2016 Ryan H. Lewis <me@ryanlewis.net> - 0.86-1
- Initial RPM
