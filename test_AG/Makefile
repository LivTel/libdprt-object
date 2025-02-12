# $Header: /space/home/eng/cjm/cvs/libdprt-object/test_AG/Makefile,v 1.1 2025-02-12 11:18:09 cjm Exp $

include ../../../Makefile.common
include ../../Makefile.common
include ../Makefile.common

INCDIR 		= $(LIBDPRT_OBJECT_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_OBJECT_DOC_HOME)/test
DOCFLAGS 	= -static
BINDIR		= $(LIBDPRT_OBJECT_BIN_HOME)/test/${HOSTTYPE}

CFLAGS 		= -g -I$(INCDIR) -I$(CFITSIOINCDIR)

SRCS 		= object_test_AG.c
OBJS 		= $(SRCS:%.c=${BINDIR}/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)

top: ${BINDIR}/object_test_AG docs

static: ${BINDIR}/object_test_AG_static docs

${BINDIR}/object_test_AG: ${BINDIR}/object_test_AG.o $(LT_LIB_HOME)/libdprt_object.so
	$(CC) -o $@ ${BINDIR}/object_test_AG.o -L/home/dev/bin/lib/i386-linux -ldprt_object -lcfitsio $(TIMELIB) -lm -lc

${BINDIR}/object_test_AG_static: ${BINDIR}/object_test_AG.o $(LT_LIB_HOME)/libdprt_object.a
	$(CC) -static -o $@ ${BINDIR}/object_test_AG.o -L/home/dev/bin/lib/i386-linux -ldprt_object -lcfitsio $(TIMELIB) -lm -lc

${BINDIR}/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

depend:
	makedepend -p$(BINDIR)/ -- $(CFLAGS) -- $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) ${BINDIR}/object_test_AG ${BINDIR}/object_test_AG_static $(OBJS) $(TIDY_OPTIONS)

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

/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/stdio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/features.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/sys/cdefs.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/wordsize.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/gnu/stubs.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/gnu/stubs-32.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/types.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/typesizes.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/libio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/_G_config.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/wchar.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/wchar.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/gconv.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/stdio_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/sys_errlist.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/stdlib.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/sys/types.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/time.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/endian.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/endian.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/sys/select.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/select.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/sigset.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/time.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/sys/sysmacros.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/pthreadtypes.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/alloca.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/string.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/math.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/huge_val.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/mathdef.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/mathcalls.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /home/dev/src/cfitsio3200/include/fitsio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/posix1_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/local_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/linux/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /usr/include/bits/posix2_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /home/dev/src/cfitsio3200/include/longnam.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test_AG.o: /home/dev/src/libdprt/object/include/object.h
