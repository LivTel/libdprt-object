# $Header: /space/home/eng/cjm/cvs/libdprt-object/c/Makefile,v 1.7 2009-01-30 15:19:48 cjm Exp $

include	../../../Makefile.common
include	../../Makefile.common
include	../Makefile.common

LIBNAME		= libdprt_object
INCDIR 		= $(LIBDPRT_OBJECT_SRC_HOME)/include
DOCSDIR 	= $(LIBDPRT_OBJECT_DOC_HOME)/cdocs
BINDIR		= $(LIBDPRT_OBJECT_BIN_HOME)/c/$(HOSTTYPE)

# log_udp library (log_udp.h is included for verbosity settings)
LOG_UDP_HOME	=	log_udp
LOG_UDP_SRC_HOME=	$(LT_SRC_HOME)/$(LOG_UDP_HOME)
LOG_UDP_CFLAGS	= -I$(LOG_UDP_SRC_HOME)/include

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
			$(LOGGINGCFLAGS) $(MEMORYCFLAGS) $(LOG_UDP_CFLAGS)
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
# Revision 1.6  2008/10/07 12:58:52  eng
# Changed DLOGGING to 6 again on CJM's instructions (jmm).
#
# Revision 1.5  2008/09/30 19:49:17  eng
# Switched off debug. Set LOGGING level to 0.
#
# Revision 1.4  2007/05/17 17:38:19  cjm
# static changes.
#
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

/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/stdarg.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/iso/stdarg_iso.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/va_list.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/stdio.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/iso/stdio_iso.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/feature_tests.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/isa_defs.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/stdio_tag.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/stdio_impl.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/stdlib.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/iso/stdlib_iso.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/string.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/iso/string_iso.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/math.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/iso/math_iso.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/time.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/types.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/machtypes.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/int_types.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/iso/time_iso.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/sys/time_impl.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /space/home/dev/src/libdprt/object/include/object.h
/space/home/dev/bin/libdprt/object/c/sun4/object.o: /usr/include/float.h
