include $(ROOTSYS)/etc/Makefile.arch
DEBUG=-g -DDEBUG
INCDIR=.

#------------------------------------------------------------------------------
PROJ = main
OBJS= $(PROJ).$(ObjSuf) singleScatterAna.$(ObjSuf)
#OBJS= $(PROJ).$(ObjSuf) neutron.$(ObjSuf)
#------------------------------------------------------------------------------
INCFLAGS = $(DEBUG) -I$(INCDIR) $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
LIBS  =  $(DEBUG) $(shell root-config --ldflags) $(ROOTLIBS) -lTreePlayer
#------------------------------------------------------------------------------
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) 
.PHONY:     $(PROJ)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(INCFLAGS) -c $< -o $@

$(PROJ):  $(OBJS)
	$(LD) $(LIBS) $^ -o $@
	@echo "$@ done"

clean:
	@rm -f $(OBJS) $(PROJ)
