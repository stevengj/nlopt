---
# NLopt on Windows
---

NLopt on Windows
----------------

NLopt works fine on Microsoft Windows computers. To simplify installation, we provide a precompiled 32-bit and 64-bit Windows [DLLs](https://en.wikipedia.org/wiki/Dynamic-link_library), built with [MinGW](https://en.wikipedia.org/wiki/MinGW):

-   [nlopt-2.4.2-dll32.zip](http://ab-initio.mit.edu/nlopt/nlopt-2.4.2-dll32.zip) (32-bit)
-   [nlopt-2.4.2-dll64.zip](http://ab-initio.mit.edu/nlopt/nlopt-2.4.2-dll64.zip) (64-bit)

Be sure to read the `README-WINDOWS` file included in this [zip](https://en.wikipedia.org/wiki/ZIP_(file_format)) archive for how to build an [import library](http://msdn.microsoft.com/en-us/library/0b9xe492.aspx) for the DLL before using, by running the `lib` `/def:libnlopt-0.def` command (which comes with Microsoft compilers). If you are using GNU compilers (MinGW), run `dlltool` `--input-def` `libnlopt-0.def` `--dllname` `libnlopt-0.dll` `--output-lib` `libnlopt-0.lib` (dlltool comes with MinGW).

(For source code, download the [main `.tar.gz` package](NLopt#Download_and_installation.md).)

Alternatively, you can use the following files, provided by Benoit Scherrer (`benoitscherrer` ατ `gmail.com`), to compile NLopt from source on Windows (with the Microsoft compiler) using [CMake](https://en.wikipedia.org/wiki/CMake):

-   [CMakeLists.txt](http://ab-initio.mit.edu/nlopt/CMakeLists.txt) and [config.cmake.h.in](http://ab-initio.mit.edu/nlopt/config.cmake.h.in)

Unofficial Python binaries for Windows are available from Christoph Gohike:

-   [Python binaries for 64-bit Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs/#nlopt)

### Octave plugin

To build the NLopt plugin for [GNU Octave](https://en.wikipedia.org/wiki/GNU_Octave) (a free Matlab clone, which uses the [same NLopt interface as in Matlab](NLopt_Matlab_Reference.md)), you will need the following additional steps. (See [Octave for Windows](http://wiki.octave.org/wiki.pl?OctaveForWindows) on the Octave web page to download Octave.)

1.  First, download *both* the `.zip` file above and the main `.tar.gz` package of NLopt. Use `dlltool` (which comes with Octave in the *octave*`\mingw23\bin` directory) to create the `.lib` import library as explained above.
2.  Unpack the NLopt .tar.gz file, and copy the `nlopt_optimize-oct.cc` from the *nlopt*`\octave` directory, along with downloading [nlopt_optimize_usage.h](http://jdj.mit.edu/~stevenj/nlopt_optimize_usage.h), and put them together with your `.lib` and `.dll` files (from the `.zip`).
3.  Compile the Octave plugin (`.oct` file) with `mkoctfile` `-lnlopt-0` `--output` `nlopt_optimize.oct` `nlopt_optimize-oct.cc` (`mkoctfile` is a program included with Octave).
4.  Finally, move `libnlopt-0.dll` to the *octave*`\bin` directory so that Octave can find it.

[Category:NLopt](index.md)
