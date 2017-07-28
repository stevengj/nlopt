---
# NLopt Installation
---

The installation of NLopt is fairly standard and straightforward, at least on Unix-like systems (GNU/Linux is fine). It doesn't require any particular packages to be installed except for a C compiler, although you need to have [Octave](https://en.wikipedia.org/wiki/GNU_Octave) and/or Matlab installed if you want to install the Octave and/or Matlab plugins, respectively.

In particular, NLopt uses the standard [Autoconf](https://en.wikipedia.org/wiki/GNU_Autoconf) `configure` script, which means that you compile it via:

```
./configure
make
```


in the `nlopt` directory. Then, you would switch to be the `root` user, or use the `sudo` command, to install the NLopt libraries and header files via:

```
make install
```


By default, this installs the NLopt static library (`libnlopt.a`) in `/usr/local/lib` and the NLopt header file (`nlopt.h`) in `/usr/local/include`, as well manual pages and a few other files.

In the following, we describe a few details of this installation process, including how to change the installation location.

Changing the installation directory
-----------------------------------

You may wish to install NLopt in a directory other than `/usr/local`, especially if you do not have administrator access to your machine. You can do this using the `--prefix` argument to the `configure` script.

For example, suppose that you want to install into the `install` subdirectory of your home directory (`$HOME`). You would do:

```
./configure --prefix=$HOME/install
make
make install
```


This will create the directories `$HOME/install/lib` etcetera and install NLopt into them. However, now when you compile code using NLopt, you will need to tell the compiler where to find the NLopt header files (using `-I`) and libraries (using `-L`) with something like:

```
cc -I$HOME/install/include myprogram.c -L$HOME/install/lib -lnlopt -lm -o myprogram
```


See also below for how to change the installation directories for Octave, Matlab, and Guile plugins, if you are installing those.

Note also that the `--prefix` flag will change the location where the Python plugins are installed, so you may need to change the [Python module search path](http://docs.python.org/tutorial/modules.html#the-module-search-path) via the `PYTHONPATH` environment variable.

Shared libraries
----------------

By default, NLopt compiles as a static library. This means that each program you link to NLopt will make a separate copy of the library, which wastes a little disk space. The alternative is to compile NLopt as a shared library (also called a dynamic-link library). While more efficient in terms of disk space etcetera, shared libraries require a bit more care to install properly, which is why we don't install them by default.

Compiling NLopt as a shared library is easy. Just add `--enable-shared` to the `configure` flags, as in:

```
./configure --enable-shared
```


Then you run `make` and `make` `install` as usual.

However, at this point you need to tell the operating system where to find the shared library, so that the runtime linker works properly. There are at least two ways to do this. First, you can use the `LD_LIBRARY_PATH` environment variable. For example, if you installed into the `/foo/bar` directory, so that the library is in `/foo/bar/lib`, then you would do

```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/foo/bar/lib
```


in the [bash](https://en.wikipedia.org/wiki/Bash) shell, or

```
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/foo/bar/lib
```


in [csh](https://en.wikipedia.org/wiki/csh) or [w:tcsh](https://en.wikipedia.org/wiki/tcsh).

Alternatively, in GNU/Linux systems, you can add the library directory to the system-wide file `/etc/ld.so.conf` and then (as root) run `/sbin/ldconfig`. [Category:NLopt](index.md)

Octave and Matlab plugins
-------------------------

When you compile NLopt using the above commands, it will automatically compile plugins for both Matlab and GNU Octave (a free Matlab clone) if the latter programs are installed. On most current systems, Matlab and Octave plugins require NLopt to be compiled as a shared library (see above).

### Matlab

In particular, for Matlab plugins to be installed, you should to have the Matlab `mex` compiler command in your [Unix PATH](http://kb.iu.edu/data/acar.html). Alternatively, you can specify the explicit path to the `mex` by passing a `MEX` variable to `configure`, via:

`./configure MEX=`*`/path/to/mex`*

Some versions of Matlab also require that you compile NLopt as a shared library in order to produce a Matlab plugin; see below.

By default, the Matlab plugins (along with help files and other `.m` files) are installed into the `MATLABPATH` printed out by `matlab` `-n` (which gives some directory within your Matlab installation directory), so that they will be available to all Matlab users. (This requires `matlab` to be in your `PATH` too; alternatively, you can pass `MATLAB=/path/to/matlab` to `configure`.) You can override this default (e.g. if you don't have administrator access on your machine) by passing a `MEX_INSTALL_DIR` to `configure`, via (in addition to other `configure` arguments):

`./configure MEX_INSTALL_DIR=`*`dir`*

to install the Matlab plugins in directory *dir*. In this case, however, when you run Matlab you will either need to run in the *dir* directory or explicitly add *dir* to your Matlab path (see the Matlab `path` command).

### Octave

For the Octave plugins to be installed, you need to have the Octave `mkoctfile` program in your PATH. `mkoctfile` is Octave's equivalent of `mex`. If you are using a GNU/Linux system, and you installed Octave using one of the precompiled packages for your distribution, then you probably need to install a *separate package* to get `mkoctfile`. For example, on Debian you need to install the `octave-headers` package, and on Redhat you need the `octave-devel` package.

By default, the compiled Octave plugins (`.oct` files) are installed into the systemwide `site/oct` directory (usually something like `/usr/lib/octave/2.1.73/site/oct/i486-pc-linux-gnu`), and the .m script files are installed into the systemwide `site/m` directory (usually something like `/usr/share/octave/2.1.73/site/m/`). You can change these defaults by passing `OCT_INSTALL_DIR` and `M_INSTALL_DIR`, respectively, to the configure script, via:

`./configure OCT_INSTALL_DIR=`*`octdir`*` M_INSTALL_DIR=`*`mdir`*

(If you only pass `OCT_INSTALL_DIR`, the default `M_INSTALL_DIR=$OCT_INSTALL_DIR`.) In this case, however you will either need to run Octave in the directory where these files are installed or explicitly add those directories to the Octave path (see the Octave `path` command).

Python plugins
--------------

If [Python](https://en.wikipedia.org/wiki/Python_(programming_language)) is installed on your machine, and you configured NLopt as a shared library (see above), then NLopt will automatically compile and install a Python `nlopt` module. You also need [NumPy](https://en.wikipedia.org/wiki/NumPy) to be installed, as NLopt's Python interface uses NumPy array types.

To specify a particular version or location of Python, use the `PYTHON` variable to set the name of the `python` executable:

`./configure PYTHON=`*`python`*

GNU Guile plugins
-----------------

If [Guile](https://en.wikipedia.org/wiki/GNU_Guile) is installed on your machine, and you configured NLopt as a shared library (see above), then a Guile `nlopt` module will automatically be compiled and installed.

Note that many GNU/Linux distributions come with only the Guile program and shared libraries pre-installed; to compile the NLopt plugin you will also need the Guile programming header files, which are usually in a `guile-dev` or `guile-devel` package that you must install separately.

If you want to specify a particular version or a nonstandard location of Guile, you should use the `GUILE_CONFIG` and `GUILE` variables to specify the locations of the `guile-config` and `guile` programs:

`./configure GUILE=`*`guile`*` GUILE_CONFIG=`*`guile-config`*

(The `configure` script uses these programs to determine the compiler flags and installation directories for Guile plugins.)

By default, the Guile plugin is installed into the system-wide site directory (the value of the `(%site-dir)` function in Guile), typically something like `/usr/share/guile/site`. If you want to install somewhere else (e.g. if you do not have administrator access), you can specify a different directory by setting `GUILE_INSTALL_DIR` on the `configure` command line:

`./configure GUILE_INSTALL_DIR=`*`dir`*

Note, however, that if you do this then Guile may not know where to load the `nlopt` module from. You can [update the Guile load path](http://www.gnu.org/software/guile/manual/html_node/Build-Config.html) by changing the `%load-path` variable in Guile or using the `GUILE_LOAD_PATH` environment variable.

NLopt with C++ algorithms
-------------------------

NLopt, as-is, is callable from C, C++, and Fortran, with optional Matlab and GNU Octave plugins (and even installs an `nlopt.hpp` C++ header file to allow you to call it in a more C++ style). By default, it includes only subroutines written in C (or written in Fortran and converted to C), to simplify linking. If you configure with:

```
./configure --with-cxx
```


however, it will also include algorithms implemented in C++ (currently, just the StoGO algorithm), and the resulting library will be called `libnlopt_cxx` and is linked with `-lnlopt_cxx`.

The `libnlopt_cxx` has the *same* interface as the ordinary NLopt library, and can *still* be called from ordinary C and Fortran programs. However, to use it you must also *link* with the C++ standard libraries. The easiest way to do this is to link with the C++ linker: compile your source files into `.o` object files, and then call the C++ compiler to link these `.o` files with `-lnlopt_cxx` into your executable program.

It is because this linking process is somewhat annoying, and it only adds a single more algorithm (StoGO) to NLopt, that by default we omit StoGO to create a library that does not require the C++ standard libraries to link.
