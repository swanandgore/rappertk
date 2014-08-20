DEPTH = ../
SUBDIRS =

CSOURCES = RanGen verbosity
ISOURCES = misc
SWIGLIB = misc

-include $(DEPTH)/makedefs
# DO NOT DELETE

RanGen.o: RanGen.h
verbosity.o: verbosity.h

misc_wrap.o: verbosity.h RanGen.h
