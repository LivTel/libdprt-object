# $Header: /space/home/eng/cjm/cvs/libdprt-object/test/Makefile,v 1.3 2012-11-12 16:22:38 cjm Exp $

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
	$(CC) -o $@ ${BINDIR}/object_test.o -L/home/dev/bin/lib/i386-linux -ldprt_object -lcfitsio $(TIMELIB) -lm -lc

${BINDIR}/object_test_static: ${BINDIR}/object_test.o $(LT_LIB_HOME)/libdprt_object.a
	$(CC) -static -o $@ ${BINDIR}/object_test.o -L/home/dev/bin/lib/i386-linux -ldprt_object -lcfitsio $(TIMELIB) -lm -lc

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
# Revision 1.2  2004/08/06 14:15:21  cjm
# New makedepend path option.
#
# Revision 1.1  2004/01/26 15:20:23  cjm
# Initial revision
#
#
# DO NOT DELETE

/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/stdio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/features.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/sys/cdefs.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/wordsize.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/gnu/stubs.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/gnu/stubs-32.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/lib/gcc/i386-redhat-linux/4.1.1/include/stddef.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/types.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/typesizes.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/libio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/_G_config.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/wchar.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/wchar.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/gconv.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/lib/gcc/i386-redhat-linux/4.1.1/include/stdarg.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/stdio_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/sys_errlist.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/stdlib.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/sys/types.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/time.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/endian.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/endian.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/sys/select.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/select.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/sigset.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/time.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/sys/sysmacros.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/pthreadtypes.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/alloca.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/strings.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/math.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/huge_val.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/mathdef.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/mathcalls.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /home/dev/src/cfitsio3310/include/fitsio.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/lib/gcc/i386-redhat-linux/4.1.1/include/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/lib/gcc/i386-redhat-linux/4.1.1/include/syslimits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/posix1_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/local_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/linux/limits.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /usr/include/bits/posix2_lim.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /home/dev/src/cfitsio3310/include/longnam.h
/home/dev/bin/libdprt/object/test/i386-linux/object_test.o: /home/dev/src/libdprt/object/include/object.h
