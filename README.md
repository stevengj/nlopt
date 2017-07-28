[![Latest Docs](https://readthedocs.org/projects/pip/badge/?version=latest)](http://nlopt.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/stevengj/nlopt.svg?branch=master)](https://travis-ci.org/stevengj/nlopt)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/stevengj/nlopt?branch=master&svg=true)](https://ci.appveyor.com/project/StevenGJohnson/nlopt)

NLopt is a library for nonlinear local and global optimization, for
functions with and without gradient information.  It is designed as
a simple, unified interface and packaging of several free/open-source
nonlinear optimization libraries.

The latest release and a complete manual may be found at the NLopt
home page: http://ab-initio.mit.edu/nlopt

NLopt is compiled and installed with the [CMake][1] build system
(see `CMakeLists.txt` file for available options):

    git clone git://github.com/stevengj/nlopt
    cd nlopt
    cmake .
    make
    sudo make install

(To build the latest development sources from git, you will need [SWIG][2]
to generate the Python and Guile bindings.)

Once it is installed, `#include <nlopt.h>` in your C/C++ programs and
link it with `-lnlopt -lm`.  You may need to use a C++ compiler to link
in order to include the C++ libraries (which are used internally by NLopt,
even though it exports a C API).

The minimization function, `nlopt_minimize`, is described in the [manpage][3]
(`api/nlopt_minimize.3`, which is installed by `make install`).
See also the manual on our web page.

There are also interfaces for Fortran, Python, MATLAB, GNU Octave, OCaml,
GNU Guile, GNU R, Lua, and Julia.  Interfaces for other languages may
be added in the future.

[1]: https://cmake.org/
[2]: http://www.swig.org/
[3]: https://en.wikipedia.org/wiki/Man_page
