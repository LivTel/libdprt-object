# $Header: /space/home/eng/cjm/cvs/libdprt-object/c/Makefile,v 1.4 2007-05-17 17:38:19 cjm Exp $

include	../../../Makefile.common
include	../../Makefile.common
include	../Makefile.common

LIBNAME		= libdprt_object
INCDIR 		= $(LIBDPRT_OBJECT_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_OBJECT_DOC_HOME)/cdocs
BINDIR		= $(LIBDPRT_OBJECT_BIN_HOME)/c/$(HOSTTYPE)

# Controls how much logging code to put into the library.
# 0 no logging
# 1 general logging
# 6 FWHM and Object logging
# 10 per pixel logging
# 11 some extra memory logging
# 6 seems about right.
LOGGINGCFLAGS	= -DLOGGING=6
# If you're insane, disable NULL pointer checks
MEMORYCFLAGS	= -DMEMORYCHECK
CFLAGS 		= -g $(CCHECKFLAG) -I$(INCDIR) -I$(CFITSIOINCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR) -L$(LT_LIB_HOME) \
			$(LOGGINGCFLAGS) $(MEMORYCFLAGS)
LINTFLAGS 	= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 	= -static
SRCS 		= object.c
HEADERS		= $(SRCS:%.c=%.h)
OBJS		= $(SRCS:%.c=$(BINDIR)/%.o)
DOCS 		= $(SRCS:%.c=$(DOCSDIR)/%.html)
LIBS		=

top: shared docs

shared: $(LT_LIB_HOME)/$(LIBNAME).so

$(LT_LIB_HOME)/$(LIBNAME).so: $(OBJS)
	$(CC) $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB) $(LIBS) $(SOCKETLIB)

$(BINDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

static: $(LT_LIB_HOME)/$(LIBNAME).a

$(LT_LIB_HOME)/$(LIBNAME).a: $(OBJS)
	ar rcv $@ $?

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

# dont checkout DpRtLibrary.h - it is a machine built header
checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)
	cd $(INCDIR); $(CO) $(CO_OPTIONS) $(HEADERS);

# dont checkin DpRtLibrary.h - it is a machine built header
checkin:
	-$(CI) $(CI_OPTIONS) $(SRCS)
	-(cd $(INCDIR); $(CI) $(CI_OPTIONS) $(HEADERS);)

staticdepend:
	makedepend -p$(BINDIR)/ -- $(CFLAGS) -- $(SRCS)

depend:
	makedepend -p$(BINDIR)/ -- $(CFLAGS) -- $(SRCS)

lint:
	$(LINT)	$(LINTFLAGS) $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) $(OBJS) $(LT_LIB_HOME)/$(LIBNAME).so $(LT_LIB_HOME)/$(LIBNAME).a 
	-$(RM) $(RM_OPTIONS) $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

backup: tidy
	-$(RM) $(RM_OPTIONS) $(OBJS)

#
# $Log: not supported by cvs2svn $
# Revision 1.3  2004/08/06 14:14:29  cjm
# New makedepend path option.
#
# Revision 1.2  2004/01/29 12:28:00  cjm
# Reduced and commented LOGGING level, performance enhancement.
# Added MEMORYCHECK CFLAGS for the insane.
#
# Revision 1.1  2004/01/26 15:20:11  cjm
# Initial revision
#
#
# DO NOT DELETE

/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/include/stdarg.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/stdio.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/features.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/sys/cdefs.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/gnu/stubs.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/lib/gcc-lib/i386-redhat-linux/3.2.2/include/stddef.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/types.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/wordsize.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/typesizes.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/libio.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/_G_config.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/wchar.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/wchar.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/gconv.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/stdio_lim.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/sys_errlist.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/stdlib.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/string.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/math.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/huge_val.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/mathdef.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/mathcalls.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/time.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /usr/include/bits/time.h
/home/dev/bin/libdprt/object/c/i386-linux/object.o: /home/dev/src/libdprt/object/include/object.h
