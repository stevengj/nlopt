#!/bin/sh
set -ev

rm -rf mingw64
make distclean || true

echo "COMPILING..."

./configure --prefix=`pwd`/mingw64 --host=x86_64-w64-mingw32 --enable-shared --disable-static --without-matlab --without-octave --without-python --without-guile --without-threadlocal && make -j4 && make install

echo "POST-PROCESSING..."

cd mingw64/bin
for dll in *.dll; do
    def=`basename $dll .dll`.def
    echo "LIBRARY $dll" > $def
    echo EXPORTS >> $def
    x86_64-w64-mingw32-nm $dll | grep ' T ' | sed 's/.* T //' | egrep 'nlopt|nlo_' >> $def
done
cd ../..

perl -pi -e 's,^ * #define NLOPT_DLL,*/\n#define NLOPT_DLL\n/*,' mingw64/include/nlopt.h

cat > README-WINDOWS <<EOF
This .zip archive contains DLL libraries and the associated header (.h)
and module-definition (.def) files of NLopt compiled for Win64.

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
x86_64-w64-mingw32-gcc --version |head -1 >> README-WINDOWS

# grep -v "nlopt-util.h" octave/nlopt_minimize_constrained-mex.c > mingw64/nlopt_minimize_constrained.c

nlopt_vers=`grep PACKAGE_VERSION config.h |cut -d" " -f3 |tr -d \"`

mkdir mingw64/matlab
cd octave
cp `grep 'MFILES =' Makefile.am | cut -d= -f2` ../mingw64/matlab
cp `grep 'm_DATA =' Makefile.am | cut -d\) -f2` ../mingw64/matlab
cp nlopt_optimize-mex.c ../mingw64/matlab/nlopt_optimize.c
cd ..

mkdir mingw64/python
cp swig/nlopt.py swig/nlopt-python.cpp mingw64/python
cat > mingw64/python/setup.py <<EOF
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
zip=nlopt-${nlopt_vers}-dll64.zip
rm -f $zip
zip -vj $zip mingw64/bin/*.dll mingw64/bin/*.exe
zip -vjgl $zip mingw64/bin/*.def mingw64/include/* mingw64/python/* README COPYING COPYRIGHT NEWS README-WINDOWS

echo "PACKAGING $zip..."

cd mingw64
zip -vgl ../$zip matlab/*
cd ..
