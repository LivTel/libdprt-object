# $Header: /space/home/eng/cjm/cvs/libdprt-object/c/Makefile,v 1.2 2004-01-29 12:28:00 cjm Exp $

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
	makedepend -- $(CFLAGS) -- $(SRCS)

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

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
# Revision 1.1  2004/01/26 15:20:11  cjm
# Initial revision
#
#
