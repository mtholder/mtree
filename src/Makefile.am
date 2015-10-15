AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

bin_PROGRAMS = mtree

noinst_HEADERS = \
	mt_char_model.h \
	mt_data.h \
	mt_log.h \
	mt_ini_options.h \
	mt_optimize_branches.h \
	mt_tree.h \
	mt_tree_traversal.h \
	INIReader.h \
	pattern_class.h \
	parsimony.h \
	ini.h \
	search.h

mtree_SOURCES = \
	likelihood.cpp \
	optimize_branches.cpp \
	main.cpp \
	ini.c \
	INIReader.cpp \
	pattern_class.cpp \
	parsimony.cpp \
	mt_opt_model.cpp \
	search.cpp

mtree_LDADD =  $(LDADD) @MTREE_LIBS@
