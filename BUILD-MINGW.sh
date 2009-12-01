./configure --prefix=`pwd`/mingw --host=i586-mingw32msvc --enable-shared --disable-static --without-matlab --without-octave && make -j4 && make install

cd mingw/bin
for dll in *.dll; do
    def=`basename $dll .dll`.def
    echo "LIBRARY $dll" > $def
    echo EXPORTS >> $def
    i586-mingw32msvc-nm $dll | grep ' T _' | sed 's/.* T _//' | grep nlopt >> $def
done
cd ../..

perl -pi -e 's,^ * #define NLOPT_DLL,*/\n#define NLOPT_DLL\n/*,' mingw/include/nlopt.h

cat > README-WINDOWS <<EOF
This .zip archive contains DLL libraries and the associated header (.h)
and module-definition (.def) files of NLopt compiled for Win32.

In order to link to this .dll files from Visual C++, you need to
create a .lib "import libraries" for it, and can do so with the "lib"
command that comes with VC++.  In particular, run:
     lib /def:libnlopt-0.def

To compile the Matlab plugin, use the Matlab "mex" compiler on the file
nlopt_minimize_constrained.c (being sure to link to the libnlopt DLL).

They were compiled by the GNU C compiler for MinGW, specifically:
EOF
i586-mingw32msvc-gcc --version |head -1 >> README-WINDOWS

# grep -v "nlopt-util.h" octave/nlopt_minimize_constrained-mex.c > mingw/nlopt_minimize_constrained.c

nlopt_vers=`grep PACKAGE_VERSION config.h |cut -d" " -f3 |tr -d \"`
zip=nlopt-${nlopt_vers}-dll.zip
rm -f $zip
zip -vj $zip mingw/bin/*.dll mingw/bin/*.exe octave/*.m mingw/nlopt_minimize_constrained.c
zip -vjgl $zip mingw/bin/*.def mingw/include/* README COPYING COPYRIGHT NEWS README-WINDOWS
