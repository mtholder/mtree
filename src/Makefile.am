AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@
 
bin_PROGRAMS = mtree

noinst_HEADERS = \
	mt_char_model.h \
	mt_data.h \
	mt_log.h \
	mt_optimize_branches.h \
	mt_tree.h \
	mt_tree_traversal.h

mtree_SOURCES = \
	likelihood.cpp \
	optimize_branches.cpp \
	main.cpp

mtree_LDADD =  $(LDADD) @MTREE_LIBS@
