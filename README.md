My modification of cscope for win32.

### Build with MinGW

    1. Install MinGw and MSYS, make sure bison (or flex and yacc) installed;
    2. Install mingw-libpdcurses;
    3. Copy 3rd include dir and lib dir into MinGW include and lib paths;
    4. in MSYS shell, run ./configure.msys and ./make.msys
