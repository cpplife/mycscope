My modification of cscope for win32.

### What done:
    1. No external sort tool; merged sort function with code;
    2. Speedup the text search with multi threads;
    3. Ignored struct declaration when search symbol.

### Build with MinGW

    1. Install MinGw and MSYS, make sure bison (or flex and yacc) installed;
    2. Install mingw-libpdcurses;
    3. Copy 3rd include dir and lib dir into MinGW include and lib paths;
    4. in MSYS shell, run ./configure.msys and ./make.msys
