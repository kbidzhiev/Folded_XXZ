#################################################################
#################################################################
#################################################################


include ./this_dir.mk
include ./options.mk


#Define Flags ----------

TENSOR_HEADERS=$(PREFIX)/itensor/core.h
CCFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(CPPFLAGS) $(OPTIMIZATIONS)
CCGFLAGS= -I. $(ITENSOR_INCLUDEFLAGS) $(DEBUGFLAGS)
LIBFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS=-L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)



#Rules ------------------

%.o: %.cc $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(ITENSOR_GLIBS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------



3siteHam.exe: 3siteHam.cc observables.cc  $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	$(CCCOM) $(CCFLAGS) 3siteHam.cc observables.cc -o $@ $(LIBFLAGS)


all:  3siteHam.exe
clean:
	rm *.exe *.o 

