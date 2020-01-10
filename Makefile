ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTCINT=$(ROOTSYS)/bin/rootcint

CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(ROOTSYS)/htmldoc -I. -O2 -g -Wall -Wshadow -W -Woverloaded-virtual -std=c++11 -fPIC $(ROOTCFLAGS)
LDFLAGS=$(ROOTLDFLAGS) -L. -Wl,-rpath .

ROOTLIBS     := -lXMLParser $(ROOTLIBS)

TUNFOLDSOURCE := TUnfoldV17.cxx TUnfoldSysV17.cxx TUnfoldDensityV17.cxx TUnfoldBinningV17.cxx TUnfoldBinningXMLV17.cxx

TUNFOLDDICT := $(subst .cxx,Dict.cxx,$(TUNFOLDSOURCE))

%Dict.cxx: %.cxx
	rm -f $@ $(@:.cxx=.h)
	$(ROOTCINT) $@ -c $<

exercise13prog: exercise13prog.C $(TUNFOLDDICT)
	$(CXX) $(CXXFLAGS) $< -o  $@ $(TUNFOLDDICT) $(LDFLAGS) $(ROOTLIBS)
