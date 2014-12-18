AM_CPPFLAGS = @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@
 
bin_PROGRAMS = mtree

noinst_HEADERS = \
	mt_char_model.h \
	mt_log.h \
	mt_tree.h \
	mt_tree_traversal.h

mtree_SOURCES = \
	likelihood.cpp \
	main.cpp

mtree_LDADD =  $(LDADD) @MTREE_LIBS@
