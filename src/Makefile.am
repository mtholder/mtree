AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

bin_PROGRAMS = mtree

noinst_HEADERS = \
	mt_char_model.h \
	mt_data.h \
	mt_log.h \
	mt_ini_options.h \
	mt_optimize_branches.h \
	mt_part_model_info.h \
	mt_tree.h \
	mt_tree_traversal.h \
	INIReader.h \
	ini.h \
  search.h

mtree_SOURCES = \
	likelihood.cpp \
	mt_optimize_branches.cpp \
	mt_part_data.cpp \
	mt_part_model_info.cpp \
	mt_prune_op_struct.cpp \
	main.cpp \
	ini.c \
	INIReader.cpp \
  search.cpp
mtree_LDADD =  $(LDADD) @MTREE_LIBS@
