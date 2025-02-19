# $Header: /space/home/eng/cjm/cvs/libdprt-object/test/Makefile,v 1.4 2025-02-12 11:38:54 cjm Exp $

include ../../../Makefile.common
include ../../Makefile.common
include ../Makefile.common

INCDIR 		= $(LIBDPRT_OBJECT_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_OBJECT_DOC_HOME)/test
DOCFLAGS 	= -static
BINDIR		= $(LIBDPRT_OBJECT_BIN_HOME)/test/${HOSTTYPE}

CFLAGS 		= -g -I$(INCDIR) -I$(CFITSIOINCDIR)

SRCS 		= object_test.c
OBJS 		= $(SRCS:%.c=${BINDIR}/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)

top: ${BINDIR}/object_test docs

static: ${BINDIR}/object_test_static docs

${BINDIR}/object_test: ${BINDIR}/object_test.o $(LT_LIB_HOME)/libdprt_object.so
	$(CC) -o $@ ${BINDIR}/object_test.o -L$(LT_LIB_HOME) -ldprt_object -lcfitsio $(TIMELIB) -lpthread -lm -lc

${BINDIR}/object_test_static: ${BINDIR}/object_test.o $(LT_LIB_HOME)/libdprt_object.a
	$(CC) -static -o $@ ${BINDIR}/object_test.o -L$(LT_LIB_HOME) -ldprt_object -lcfitsio $(TIMELIB) -lpthread -lm -lc

${BINDIR}/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

depend:
	makedepend $(MAKEDEPENDFLAGS) -p$(BINDIR)/ -- $(CFLAGS) -- $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) ${BINDIR}/object_test ${BINDIR}/object_test_static $(OBJS) $(TIDY_OPTIONS)

tidy:
	-$(RM) $(RM_OPTIONS) $(TIDY_OPTIONS)

backup: tidy
#	-$(RM) $(RM_OPTIONS) 

checkin:
	-$(CI) $(CI_OPTIONS) $(SRCS)

checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)

#
# $Log: not supported by cvs2svn $
# Revision 1.3  2012/11/12 16:22:38  cjm
# makedepend flags.
#
# Revision 1.2  2004/08/06 14:15:21  cjm
# New makedepend path option.
#
# Revision 1.1  2004/01/26 15:20:23  cjm
# Initial revision
#
#
