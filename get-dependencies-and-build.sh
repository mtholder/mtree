#!/bin/sh
set -x
export MTREE_ROOT="$PWD"
# build NCL
if ! test -d ncl
then
    git clone https://github.com/mtholder/ncl.git || exit 1
fi
cd ncl || exit
autoreconf || exit
CPPLAGS="-D__extern_always_inline=inline" CC=clang CXX=clang++ ./configure --prefix=$MTREE_ROOT/installed --with-constfuncs || exit
make -j4 || exit 
make install || exit
cd ..
export LD_LIBRARY_PATH="$MTREE_ROOT/installed/lib/ncl:$LD_LIBRARY_PATH";

# build configure script
autoreconf
automake --add-missing
autoreconf || exit


if ! test -d build
then
    mkdir build || exit
fi
cd build || exit
CC=clang CXX=clang++ ../configure --prefix=$MTREE_ROOT/installed --with-ncl=$MTREE_ROOT/installed --enable-debugging=yes --enable-asserts || exit
make -j4 || exit


make check || exit
#./src/mtree ./tests/binary.nex
