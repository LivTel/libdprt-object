# $Header: /space/home/eng/cjm/cvs/libdprt-object/test_AG_loop/Makefile,v 1.1 2025-02-12 11:15:36 cjm Exp $

include ../../../Makefile.common
include ../../Makefile.common
include ../Makefile.common

INCDIR 		= $(LIBDPRT_OBJECT_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_OBJECT_DOC_HOME)/test
DOCFLAGS 	= -static
BINDIR		= $(LIBDPRT_OBJECT_BIN_HOME)/test/${HOSTTYPE}

CFLAGS 		= -g -I$(INCDIR) -I$(CFITSIOINCDIR)

SRCS 		= object_test_AG_loop.c
OBJS 		= $(SRCS:%.c=${BINDIR}/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)

top: ${BINDIR}/object_test_AG_loop docs

static: ${BINDIR}/object_test_AG_loop_static docs

${BINDIR}/object_test_AG_loop: ${BINDIR}/object_test_AG_loop.o $(LT_LIB_HOME)/libdprt_object.so
	$(CC) -o $@ ${BINDIR}/object_test_AG_loop.o -L/home/dev/bin/lib/i386-linux -ldprt_object -lcfitsio $(TIMELIB) -lm -lc

${BINDIR}/object_test_AG_loop_static: ${BINDIR}/object_test_AG_loop.o $(LT_LIB_HOME)/libdprt_object.a
	$(CC) -static -o $@ ${BINDIR}/object_test_AG_loop.o -L/home/dev/bin/lib/i386-linux -ldprt_object -lcfitsio $(TIMELIB) -lm -lc

${BINDIR}/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

depend:
	makedepend -p$(BINDIR)/ -- $(CFLAGS) -- $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) ${BINDIR}/object_test_AG_loop ${BINDIR}/object_test_AG_loop_static $(OBJS) $(TIDY_OPTIONS)

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
# Revision 1.2  2004/08/06 14:15:21  cjm
# New makedepend path option.
#
# Revision 1.1  2004/01/26 15:20:23  cjm
# Initial revision
#
#
# DO NOT DELETE

/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/stdio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/features.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/sys/cdefs.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/gnu/stubs.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/include/stddef.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/types.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/wordsize.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/typesizes.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/libio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/_G_config.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/wchar.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/wchar.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/gconv.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/include/stdarg.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/stdio_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/sys_errlist.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/stdlib.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/string.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/time.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/time.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/math.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/huge_val.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/mathdef.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/bits/mathcalls.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /home/dev/src/cfitsio3.006/include/fitsio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/include/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/include/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /home/dev/src/cfitsio3.006/include/longnam.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG_loop.o: /home/dev/src/libdprt/object/include/object.h
