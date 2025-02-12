# Makefile
# $Header: /space/home/eng/cjm/cvs/libdprt-object/Makefile,v 1.3 2025-02-12 11:10:59 cjm Exp $ 

include ../../Makefile.common
include ../Makefile.common
include Makefile.common

DIRS = c test

top:
	@for i in $(DIRS); \
	do \
		(echo making in $$i...; cd $$i; $(MAKE) ); \
	done;

static:
	@for i in $(DIRS); \
	do \
		(echo making in $$i...; cd $$i; $(MAKE) static); \
	done;

checkin:
	-@for i in $(DIRS); \
	do \
		(echo checkin in $$i...; $(MAKE) -C $$i checkin; cd $$i; $(CI) $(CI_OPTIONS) Makefile); \
	done;

checkout:
	@for i in $(DIRS); \
	do \
		(echo checkout in $$i...; cd $$i; $(CO) $(CO_OPTIONS) Makefile; $(MAKE) checkout); \
	done;

depend:
	@for i in $(DIRS); \
	do \
		(echo depend in $$i...; $(MAKE) -C $$i depend); \
	done;

clean:
	$(RM) $(RM_OPTIONS) $(TIDY_OPTIONS)
	@for i in $(DIRS); \
	do \
		(echo clean in $$i...;  cd $$i; $(MAKE) clean); \
	done;

tidy:
	$(RM) $(RM_OPTIONS) $(TIDY_OPTIONS)
	@for i in $(DIRS); \
	do \
		(echo tidy in $$i...; $(MAKE) -C $$i tidy); \
	done;

backup: tidy checkin
	@for i in $(DIRS); \
	do \
		(echo backup in $$i...; $(MAKE) -C $$i backup); \
	done;
	tar cvf $(BACKUP_DIR)/libdprt.tar .
	compress $(BACKUP_DIR)/libdprt.tar

# $Log: not supported by cvs2svn $
# Revision 1.2  2007/01/09 15:35:17  cjm
# Changed 'cd's to make -C.
#
# Revision 1.1  2004/02/04 14:08:29  cjm
# Initial revision
#
# Revision 1.2  2001/05/15 16:43:57  cjm
# New libdprt version.
#
# Revision 1.1  2001/05/15 16:42:15  cjm
# Initial revision
#
# Revision 1.1  1999/12/10 12:26:20  cjm
# Initial revision
#
