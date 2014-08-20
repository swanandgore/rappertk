DEPTH = ../
SUBDIRS =

LINCFLAGS = $(RTKROOT)/geometry/_geometry.so $(RTKROOT)/samplers/_samplers.so $(RTKROOT)/misc/_misc.so

CSOURCES = Constants \
	    Builder BasicPepBuilder DihedBuilder PeptideBuilder CBbuilder ChiBuilder PeptideBridgeBuilder \
        NAbuilder NAsuiteBuilder CNCaBuilder InformedPeptideBuilder NanchorBuilder BuilderGroup RapPepBuilder \
        ProtResBuilder OptFragPlacer SCLbuilder \
        TransRotator Rotator

ISOURCES = builders
SWIGLIB = builders

-include $(DEPTH)/makedefs
# DO NOT DELETE

BasicPepBuilder.o: BasicPepBuilder.h Builder.h Constants.h
BasicPepBuilder.o: /usr/local/xtal/rtk152-O4/samplers/PhipsiSampler.h
BasicPepBuilder.o: /usr/local/xtal/rtk152-O4/samplers/Sampler.h
BasicPepBuilder.o: /usr/local/xtal/rtk152-O4/samplers/OmegaSampler.h
BasicPepBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
BasicPepBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
BasicPepBuilder.o: /usr/local/xtal/rtk152-O4/misc/RanGen.h
Builder.o: Builder.h Constants.h
BuilderGroup.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h BuilderGroup.h
BuilderGroup.o: Builder.h Constants.h
CBbuilder.o: CBbuilder.h Builder.h
CBbuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
CBbuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h Constants.h
ChiBuilder.o: ChiBuilder.h Builder.h
ChiBuilder.o: /usr/local/xtal/rtk152-O4/samplers/BBdepChiSampler.h
ChiBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
ChiBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h Constants.h
CNCaBuilder.o: CNCaBuilder.h Builder.h Constants.h
CNCaBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
CNCaBuilder.o: /usr/local/xtal/rtk152-O4/misc/RanGen.h
Constants.o: Constants.h
DihedBuilder.o: DihedBuilder.h Builder.h
DihedBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
DihedBuilder.o: /usr/local/xtal/rtk152-O4/misc/RanGen.h
InformedPeptideBuilder.o: InformedPeptideBuilder.h BasicPepBuilder.h
InformedPeptideBuilder.o: Builder.h
InformedPeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/InformedPhipsiSampler.h
InformedPeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/PhipsiSampler.h
InformedPeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/Sampler.h
InformedPeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/OmegaSampler.h
InformedPeptideBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
NAbuilder.o: NAbuilder.h Builder.h Constants.h
NAbuilder.o: /usr/local/xtal/rtk152-O4/samplers/NAsampler.h
NAbuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
NanchorBuilder.o: /usr/local/xtal/rtk152-O4/misc/RanGen.h
NanchorBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
NanchorBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
NanchorBuilder.o: NanchorBuilder.h Builder.h Constants.h
NAsuiteBuilder.o: NAsuiteBuilder.h Builder.h
NAsuiteBuilder.o: /usr/local/xtal/rtk152-O4/samplers/NAsuiteSampler.h
NAsuiteBuilder.o: Constants.h NAbuilder.h
NAsuiteBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
NAsuiteBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
OptFragPlacer.o: OptFragPlacer.h Builder.h
OptFragPlacer.o: /usr/local/xtal/rtk152-O4/samplers/FragSampler.h
OptFragPlacer.o: /usr/local/xtal/rtk152-O4/samplers/Sampler.h
OptFragPlacer.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
PeptideBridgeBuilder.o: PeptideBridgeBuilder.h Builder.h
PeptideBridgeBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
PeptideBridgeBuilder.o: /usr/local/xtal/rtk152-O4/builders/Constants.h
PeptideBridgeBuilder.o: /usr/local/xtal/rtk152-O4/samplers/RATdata.h
PeptideBridgeBuilder.o: /usr/local/xtal/rtk152-O4/misc/RanGen.h
PeptideBuilder.o: PeptideBuilder.h BasicPepBuilder.h Builder.h Constants.h
PeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/PhipsiSampler.h
PeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/Sampler.h
PeptideBuilder.o: /usr/local/xtal/rtk152-O4/samplers/OmegaSampler.h
PeptideBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
PeptideBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
ProtResBuilder.o: ProtResBuilder.h Builder.h BasicPepBuilder.h CBbuilder.h
ProtResBuilder.o: ChiBuilder.h
ProtResBuilder.o: /usr/local/xtal/rtk152-O4/samplers/PhipsiSampler.h
ProtResBuilder.o: /usr/local/xtal/rtk152-O4/samplers/Sampler.h
ProtResBuilder.o: /usr/local/xtal/rtk152-O4/samplers/OmegaSampler.h
ProtResBuilder.o: /usr/local/xtal/rtk152-O4/samplers/BBdepChiSampler.h
ProtResBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
ProtResBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
RapPepBuilder.o: RapPepBuilder.h BasicPepBuilder.h Builder.h Constants.h
RapPepBuilder.o: /usr/local/xtal/rtk152-O4/samplers/RapSampler.h
RapPepBuilder.o: /usr/local/xtal/rtk152-O4/samplers/OmegaSampler.h
RapPepBuilder.o: /usr/local/xtal/rtk152-O4/samplers/Sampler.h
RapPepBuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
RapPepBuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
Rotator.o: Rotator.h Builder.h /usr/local/xtal/rtk152-O4/geometry/functions.h
Rotator.o: /usr/local/xtal/rtk152-O4/misc/RanGen.h
SCLbuilder.o: SCLbuilder.h Builder.h
SCLbuilder.o: /usr/local/xtal/rtk152-O4/samplers/SCLsampler.h
SCLbuilder.o: /usr/local/xtal/rtk152-O4/geometry/functions.h
SCLbuilder.o: /usr/local/xtal/rtk152-O4/misc/verbosity.h
SStweakerED.o: SStweakerED.h
TransRotator.o: TransRotator.h Builder.h
TransRotator.o: /usr/local/xtal/rtk152-O4/geometry/functions.h

builders_wrap.o: Builder.h DihedBuilder.h PeptideBuilder.h BasicPepBuilder.h
builders_wrap.o: PeptideBridgeBuilder.h RapPepBuilder.h
builders_wrap.o: InformedPeptideBuilder.h CNCaBuilder.h CBbuilder.h
builders_wrap.o: ChiBuilder.h Constants.h NAbuilder.h NAsuiteBuilder.h
builders_wrap.o: NanchorBuilder.h BuilderGroup.h ProtResBuilder.h
builders_wrap.o: TransRotator.h Rotator.h OptFragPlacer.h SCLbuilder.h
