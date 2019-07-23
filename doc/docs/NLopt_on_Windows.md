---
# NLopt on Windows
---

NLopt on Windows
----------------

NLopt works fine on Microsoft Windows computers, and you can compile it directly using the included [CMake](https://en.wikipedia.org/wiki/CMake) build scripts.

To simplify installation, there are also precompiled 32-bit and 64-bit Windows [DLLs](https://en.wikipedia.org/wiki/Dynamic-link_library) (along with binaries for many other systems) at [NLoptBuilder/releases](https://github.com/stevengj/NLoptBuilder/releases).  In particular, the Windows builds are

-  [NLopt.v2.6.1.i686-w64-mingw32.tar.gz](https://github.com/stevengj/NLoptBuilder/releases/download/v2.6.1/NLopt.v2.6.1.i686-w64-mingw32.tar.gz) (32-bit)
-  [NLopt.v2.6.1.x86_64-w64-mingw32.tar.gz](https://github.com/stevengj/NLoptBuilder/releases/download/v2.6.1/NLopt.v2.6.1.x86_64-w64-mingw32.tar.gz) (64-bit)

These `.tar.gz` files unpack (with a variety of Windows software, e.g. 7-zip) into a folder with a `bin` subdirectory that contains `libnlopt.dll`.  To link with this in your compiler, you will typically also want the [import library](https://stackoverflow.com/questions/3573475/how-does-the-import-library-work-details) for the DLL, which can be found in the `lib` subdirectory and is called `libnlopt.dll.a` (this can be used similarly to the `.lib` files you may be used to).

Unofficial Python binaries for Windows are available from Christoph Gohike:

-   [Python binaries for 64-bit Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs/#nlopt)

### NLopt with MinGW

If you want to compile NLopt on Windows:

1. Install [MSYS2](https://www.msys2.org/).
2. Run the following command:
    `pacman -S --needed base-devel mingw-w64-i686-toolchain mingw-w64-x86_64-toolchain \
                    mingw-w64-i686-cmake mingw-w64-x86_64-cmake --disable-download-timeout`
3. Download the NLopt using git.
    `git clone https://github.com/stevengj/nlopt.git`
4. Create a build folder inside the project.
5. From build folder run the following command:
    `cmake -G"MSYS Makefiles" ..`
6. Then run `make`.

### Octave plugin

To build the NLopt plugin for [GNU Octave](https://en.wikipedia.org/wiki/GNU_Octave) (a free Matlab clone, which uses the [same NLopt interface as in Matlab](NLopt_Matlab_Reference.md)), you will need the following additional steps. (See [Octave for Windows](https://wiki.octave.org/Octave_for_Microsoft_Windows) on the Octave web page to download Octave.)

1. Copy `libnlopt.dll` from `build` to `C:\Octave\Octave-X.X.X.X\mingw64\bin`.
2. Copy `libnlopt.dll.a` from `build` to `C:\Octave\Octave-X.X.X.X\mingw64\lib`.
3. Change the current folder to 'src\octave' and compile the Octave plugin (`.oct` file) with `mkoctfile -I"build/src/api" -lnlopt --output nlopt_optimize nlopt_optimize-oct.cc` (`mkoctfile` is a program included with Octave).