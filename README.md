# mtree
likelihood-based morphological phylogenetics

## Installation
The INSTALL file has the automake/autoconf boilerplate. The basics
are:

    $ autreconf # may fail on the first time
    $ automake add-missing
    $ autoreconf # should succeed
    $ CC=/usr/bin/clang CXX=/usr/bin/clang++ ./configure --prefix=$PWD/installed --with-ncl=/tmpncl/install --enable-debugging=yes --enable-asserts
    $ make
    $ make check

## Running
The two arguments should be a data file (NEXUS by default) and
    a -m flag followed by the path to an INI-style configuration
    file:

    $ ./src/mtree -mtests/jc.ini jc.nex

Will trigger a scoring of the tree in jc.nex with the data in 
    that file using an Mk model.

### Configuration settings
The section and setting names of the [INI file](http://en.wikipedia.org/wiki/INI_file) are
    case-sensitive.
The values are only case-sensitive when they refer to filepaths.

The currently supported values (using section/setting as a concise
    way to indicate the section and setting name) are:

| Section/Setting | Values | Meaning |
| --------------- |--------|---------|
| action/action   | LScore | Likelihood score with no branch length or parameter optmization |


`action/action = LScore (score the tree with no branch length or parameter optimization)
## Credits
Code written by Liam Heins and Mark T. Holder
with many routines/ideas taken from:
    Nexus Class Library (by Paul Lewis)
    Phylogenetic Likelihood Library (by Flouris et at)
    PAML (by Ziheng Yang)

Test code uses bash functions by Mitch Frazier from
    http://www.linuxjournal.com/content/floating-point-math-bash
INI file parsing code (ini.c, ini.h, INIReader.cpp, INIReader.cpp)
    is from the inih project https://code.google.com/p/inih/ which
    is released under the BSD License.

Thanks the US National Science Foundation and the
    Heidelberg Institute for Theoretical Studies.

