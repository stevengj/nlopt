#!/bin/sh
set -ev

rm -rf mingw32
make distclean || true

echo "COMPILING..."

./configure --prefix=`pwd`/mingw32 --host=i686-w64-mingw32 --enable-shared --disable-static --without-matlab --without-octave --without-python --without-guile --without-threadlocal && make -j4 && make install

echo "POST-PROCESSING..."

cd mingw32/bin
for dll in *.dll; do
    def=`basename $dll .dll`.def
    echo "LIBRARY $dll" > $def
    echo EXPORTS >> $def
    i686-w64-mingw32-nm $dll | grep ' T _' | sed 's/.* T _//' | egrep 'nlopt|nlo_' >> $def
done
cd ../..

perl -pi -e 's,^ * #define NLOPT_DLL,*/\n#define NLOPT_DLL\n/*,' mingw32/include/nlopt.h

cat > README-WINDOWS <<EOF
This .zip archive contains DLL libraries and the associated header (.h)
and module-definition (.def) files of NLopt compiled for Win32.

In order to link to this .dll files from Visual C++, you need to
create a .lib "import libraries" for it, and can do so with the "lib"
command that comes with VC++.  In particular, run:
     lib /def:libnlopt-0.def

To compile the Matlab plugin, use the Matlab "mex" compiler on the file
nlopt_optimize.c (being sure to link to the libnlopt DLL) in the matlab
subdirectory.

To build the Python plugin (assuming that you have Python and Numpy
installed), do:
   python setup.py build_ext --inplace

They were compiled by the GNU C compiler for MinGW, specifically:
EOF
i686-w64-mingw32-gcc --version |head -1 >> README-WINDOWS

# grep -v "nlopt-util.h" octave/nlopt_minimize_constrained-mex.c > mingw32/nlopt_minimize_constrained.c

nlopt_vers=`grep PACKAGE_VERSION config.h |cut -d" " -f3 |tr -d \"`

mkdir mingw32/matlab
cd octave
cp `grep 'MFILES =' Makefile.am | cut -d= -f2` ../mingw32/matlab
cp `grep 'm_DATA =' Makefile.am | cut -d\) -f2` ../mingw32/matlab
cp nlopt_optimize-mex.c ../mingw32/matlab/nlopt_optimize.c
cd ..

mkdir mingw32/python
cp swig/nlopt.py swig/nlopt-python.cpp mingw32/python
cat > mingw32/python/setup.py <<EOF
from distutils.core import setup, Extension
nlopt_module = Extension('_nlopt',
                           sources=['nlopt-python.cpp'],
                           libraries=['libnlopt-0'],
                           )
import numpy
setup (name = 'nlopt',
       version = '${nlopt_vers}',
       author      = "Steven G. Johnson",
       description = """NLopt nonlinear-optimization library""",
       ext_modules = [nlopt_module],
       py_modules = ["nlopt"],
       include_dirs = ['.', numpy.get_include()],
       )
EOF

nlopt_vers=`grep PACKAGE_VERSION config.h |cut -d" " -f3 |tr -d \"`
zip=nlopt-${nlopt_vers}-dll32.zip
rm -f $zip
zip -vj $zip mingw32/bin/*.dll mingw32/bin/*.exe
zip -vjgl $zip mingw32/bin/*.def mingw32/include/* mingw32/python/* README COPYING COPYRIGHT NEWS README-WINDOWS

echo "PACKAGING $zip..."

cd mingw32
zip -vgl ../$zip matlab/*
cd ..
