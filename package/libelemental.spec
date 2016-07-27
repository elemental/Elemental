Name: elemental
Version: 0.85
Release: 1%{?dist}
Summary: A fast generic C++ library for applied and computational topology
Group: Development/Libraries
License: BSD	
URL: http://libelemental.org/documentation/ 
Source0: http://libelemental.org/pub/releases/Elemental-%{version}.tgz 
BuildRequires: cmake
BuildRequires: metis-devel >= 5.1.0
%description
Elemental is an open-source library for distributed-memory dense and sparse-direct linear algebra and optimization which builds on top of BLAS, LAPACK, and MPI using modern C++ and additionally exposes interfaces to C and Python (with a Julia interface beginning development).

%package devel
Summary: Elemental Library Headers 
Group: Development/Libraries
Requires: %{name} = %{version}-%{release}
Provides: %{name}-static = %{version}-%{release}

%description devel
C++ Header files for the Elemental packages.


%package docs
Summary: Elemental Documentation package
Group: Development/Libraries
Requires: %{name} = %{version}-%{release}

%description docs
Documentation in HTML format for Elemental. This is the
same as the documentation on the Elemental website
(http://libelemental.org/documentation/)

%package examples
Summary: Source examples for Elemental
Requires: %{name} = %{version}-%{release}
#require devel so that if a user wants to compile an example they have all the headers.

%description examples 
This package contains example source files 

%prep
%setup -q -n ctl
#this file has not been written in v0.1 
cmake -DCMAKE_BUILD_TYPE=Release ..

%build
#make (in parallel)
make %{?_smp_mflags} CFLAGS="%{optflags}"

%install
#create output locatioins
mkdir -p %{buildroot}%{_bindir}
mkdir -p %{buildroot}%{_includedir}
mkdir -p %{buildroot}%{_docdir}
mkdir -p %{buildroot}%{_mandir}/man1
mkdir -p %{buildroot}%{_docdir}/examples

#cmake does the deed
make install

#put stuff in the right directories
cp -pr bin/* %{buildroot}%{_bindir}
cp -pr include/* %{buildroot}%{_includedir}
cp -pr doc/* %{buildroot}%{_docdir}
cp -pr man/* %{buildroot}%{_mandir}/man1
cp -pr examples/* %{buildroot}%{_docdir}/examples

%files  
%{_bindir}/*
%doc %{_mandir}/man1/*

%files devel
%{_includedir}/*

%files docs
%doc %{_docdir}/

%files examples
%doc %{_docdir}/examples/*

%changelog
* Sun Apr 27 2014 Ryan H. Lewis <me@ryanlewis.net> - 0.1-1
- Initial Package 
