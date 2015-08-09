# Makefile for CDFTOOLS_3.0

# ( make.macro is a link that points to the file macro.xxx where 
#   xxx is representative of your machine )
# !!----------------------------------------------------------------------
# !! CDFTOOLS_3.0 , MEOM 2011
# !! $Id$
# !! Copyright (c) 2010, J.-M. Molines
# !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
# !!----------------------------------------------------------------------

include make.macro

BINDIR = ./bin

VPATH = $(BINDIR)

EXEC = cdfmoy cdfmoyt cdfstd  cdfmoy_weighted cdfmoy_freq cdfvT cdfuv\
       cdfvsig cdfspeed cdfsum cdfenstat\
       cdfmoyuvwt \
       cdfeke cdfrmsssh cdfstdevw cdfstdevts cdflinreg cdfimprovechk\
       cdfstats \
       cdfbn2 cdfrichardson cdfsig0 cdfsigntr cdfsigi cdfsiginsitu cdfbottomsig cdfbotpressure cdfspice\
       cdfbottom cdfets cdfokubo-w cdfcurl cdfdiv cdflap cdfw cdfgeo-uv cdfmxl \
       cdfrhoproj  cdfzisot cdfsigintegr cdfpvor \
       cdfmhst cdfvhst cdfvtrp cdftransport cdfvFWov \
       cdfsigtrp cdfsigtrp_broken cdftempvol-full\
       cdfpsi cdfmoc  cdfmocsig cdfmean cdfdegradt cdfdegradu cdfdegradv cdfdegradw \
       cdfheatc cdfzonalmean cdfzonalmeanvT cdfhflx cdfwflx cdfbuoyflx\
       cdfmxlheatc cdfmxlsaltc cdfmxlhcsc cdfvertmean cdfvint \
       cdfpendep cdfzonalsum cdficediags cdfzonalout\
       cdfprofile  cdfwhereij cdffindij cdfweight cdfmaxmoc cdfcensus cdfzoom cdfmax cdfprobe cdfinfo \
       cdf16bit cdfvita cdfvita-geo cdfconvert cdfflxconv cdfclip cdfsstconv cdfstrconv cdfbathy cdfvar \
       cdfcsp cdfcoloc cdfmltmask cdfstatcoord  cdfpolymask cdfsmooth cdfmkmask cdfdifmask\
       cdfgradT cdfeddyscale_pass1 cdfeddyscale \
       cdfkempemekeepe cdfbci cdfbti cdfnrjcomp cdfcofdis cdfsections cdfnorth_unfold cdfovide cdfmppini\
       cdf_xtrac_brokenline cdf2levitusgrid2d  cdf2levitusgrid3d\
       cdfpsi_level cdfhdy cdfhdy3d cdffracinv  cdfmaskdmp cdfnan cdfscale cdfnamelist \
       cdfisopsi cdf2matlab cdffixtime cdfgeostrophy cdfchgrid cdfcmp

.PHONY: all help clean cleanexe install man installman

all: $(EXEC)

help:
	@echo "#-------------------------------------------------"
	@echo "# List of make targets:"
	@echo "#  all          : build cdftools binary"
	@echo "#  man          : build manual"
	@echo "#  clean        : remove building object (.o, .mod...)"
	@echo "#  cleanexe     : remove binary executable"
	@echo "#  install      : install binary in INSTALL folder"
	@echo "#  install_man  : install manual in INSTALL_MAN folder"
	@echo "#-------------------------------------------------"

## Statistical programs
cdfmoy: cdfio.o   cdfmoy.f90
	$(F90) cdfmoy.f90 -o $(BINDIR)/cdfmoy cdfio.o  modcdfnames.o $(FFLAGS)

cdfmoyt: cdfio.o   cdfmoyt.f90
	$(F90) cdfmoyt.f90 -o $(BINDIR)/cdfmoyt cdfio.o modcdfnames.o  $(FFLAGS)

cdfenstat: cdfio.o   cdfenstat.f90
	$(F90) cdfenstat.f90 -o $(BINDIR)/cdfenstat cdfio.o modcdfnames.o  $(FFLAGS)

cdfenstat2: cdfio.o   cdfenstat2.f90
	$(F90) cdfenstat2.f90 -o $(BINDIR)/cdfenstat2 cdfio.o modcdfnames.o  $(FFLAGS)

cdfmoy_freq: cdfio.o   cdfmoy_freq.f90
	$(F90) cdfmoy_freq.f90 -o $(BINDIR)/cdfmoy_freq cdfio.o  modcdfnames.o $(FFLAGS)

cdfmoyuvwt: cdfio.o   cdfmoyuvwt.f90
	$(F90) cdfmoyuvwt.f90 -o $(BINDIR)/cdfmoyuvwt cdfio.o modcdfnames.o $(FFLAGS)

cdfstd: cdfio.o  cdfstd.f90
	$(F90)  cdfstd.f90 -o $(BINDIR)/cdfstd cdfio.o modcdfnames.o $(FFLAGS)

cdfmoy_weighted: cdfio.o   cdfmoy_weighted.f90
	$(F90) cdfmoy_weighted.f90 -o $(BINDIR)/cdfmoy_weighted cdfio.o modcdfnames.o $(FFLAGS)

cdfeke: cdfio.o  cdfeke.f90
	$(F90) cdfeke.f90 -o $(BINDIR)/cdfeke cdfio.o modcdfnames.o $(FFLAGS)

cdfrmsssh: cdfio.o  cdfrmsssh.f90
	$(F90) cdfrmsssh.f90 -o $(BINDIR)/cdfrmsssh cdfio.o modcdfnames.o $(FFLAGS)

cdfstdevw: cdfio.o  cdfstdevw.f90
	$(F90) cdfstdevw.f90 -o $(BINDIR)/cdfstdevw cdfio.o modcdfnames.o $(FFLAGS)

cdfstdevts: cdfio.o  cdfstdevts.f90
	$(F90) cdfstdevts.f90 -o $(BINDIR)/cdfstdevts cdfio.o modcdfnames.o $(FFLAGS)

cdfvT: cdfio.o  modutils.o cdfvT.f90
	$(F90) cdfvT.f90 -o $(BINDIR)/cdfvT cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfuv: cdfio.o  modutils.o cdfuv.f90
	$(F90) cdfuv.f90 -o $(BINDIR)/cdfuv cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfvsig: cdfio.o eos.o modutils.o  cdfvsig.f90
	$(F90) $(OMP) cdfvsig.f90 -o $(BINDIR)/cdfvsig cdfio.o eos.o modcdfnames.o modutils.o $(FFLAGS)

cdfspeed: cdfio.o  cdfspeed.f90
	$(F90) cdfspeed.f90 -o $(BINDIR)/cdfspeed cdfio.o modcdfnames.o $(FFLAGS)

cdfspeedlog: cdfio.o  cdfspeedlog.f90
	$(F90) cdfspeedlog.f90 -o $(BINDIR)/cdfspeedlog cdfio.o modcdfnames.o $(FFLAGS)

cdfimprovechk: cdfio.o  cdfimprovechk.f90
	$(F90) cdfimprovechk.f90 -o $(BINDIR)/cdfimprovechk cdfio.o modcdfnames.o $(FFLAGS)

cdfstats: cdfio.o  cdfstats.f90
	$(F90) cdfstats.f90 -o $(BINDIR)/cdfstats cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdflinreg: cdfio.o  cdflinreg.f90
	$(F90) cdflinreg.f90 -o $(BINDIR)/cdflinreg cdfio.o modcdfnames.o $(FFLAGS)

## Derived quantities programs
cdfbn2: cdfio.o  eos.o  cdfbn2.f90
	$(F90) $(OMP) cdfbn2.f90 -o $(BINDIR)/cdfbn2 cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfrichardson: cdfio.o  eos.o  cdfrichardson.f90
	$(F90) $(OMP) cdfrichardson.f90 -o $(BINDIR)/cdfrichardson cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfsig0: cdfio.o  eos.o  cdfsig0.f90
	$(F90) $(OMP) cdfsig0.f90 -o $(BINDIR)/cdfsig0 cdfio.o eos.o modcdfnames.o  $(FFLAGS)

cdfsigntr: cdfio.o  eos.o  cdfsigntr.f90
	$(F90) $(OMP) cdfsigntr.f90 -o $(BINDIR)/cdfsigntr cdfio.o eos.o modcdfnames.o  $(FFLAGS)

cdfspice: cdfio.o  eos.o  cdfspice.f90
	$(F90) $(OMP) cdfspice.f90 -o $(BINDIR)/cdfspice cdfio.o eos.o modcdfnames.o  $(FFLAGS)

cdfsigi: cdfio.o  eos.o  cdfsigi.f90
	$(F90) $(OMP) cdfsigi.f90 -o $(BINDIR)/cdfsigi cdfio.o eos.o  modcdfnames.o $(FFLAGS)

cdfsiginsitu: cdfio.o  eos.o  cdfsiginsitu.f90
	$(F90) $(OMP) cdfsiginsitu.f90 -o $(BINDIR)/cdfsiginsitu cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfbottomsig: cdfio.o  eos.o  cdfbottomsig.f90
	$(F90) $(OMP) cdfbottomsig.f90 -o $(BINDIR)/cdfbottomsig cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfbotpressure: cdfio.o  eos.o modutils.o  cdfbotpressure.f90
	$(F90) $(OMP) cdfbotpressure.f90 -o $(BINDIR)/cdfbotpressure cdfio.o eos.o modcdfnames.o modutils.o $(FFLAGS)

cdfbottom: cdfio.o    cdfbottom.f90
	$(F90) cdfbottom.f90 -o $(BINDIR)/cdfbottom cdfio.o  modcdfnames.o $(FFLAGS)

cdfets: cdfio.o  eos.o  cdfets.f90
	$(F90) $(OMP) cdfets.f90 -o $(BINDIR)/cdfets cdfio.o eos.o modcdfnames.o  $(FFLAGS)

cdfokubo-w: cdfio.o  cdfokubo-w.f90
	$(F90) cdfokubo-w.f90 -o $(BINDIR)/cdfokubo-w cdfio.o modcdfnames.o $(FFLAGS)

cdfmsk: cdfio.o  cdfmsk.f90
	$(F90) cdfmsk.f90 -o $(BINDIR)/cdfmsk cdfio.o modcdfnames.o $(FFLAGS)

cdfmkmask: cdfio.o  cdfmkmask.f90
	$(F90) cdfmkmask.f90 -o $(BINDIR)/cdfmkmask cdfio.o modcdfnames.o $(FFLAGS)

cdfmeshmask: cdfio.o  cdfmeshmask.f90
	$(F90) cdfmeshmask.f90 -o $(BINDIR)/cdfmeshmask cdfio.o modcdfnames.o $(FFLAGS)

cdfmltmask: cdfio.o  cdfmltmask.f90
	$(F90) cdfmltmask.f90 -o $(BINDIR)/cdfmltmask cdfio.o modcdfnames.o $(FFLAGS)

cdfdifmask: cdfio.o  cdfdifmask.f90
	$(F90) cdfdifmask.f90 -o $(BINDIR)/cdfdifmask cdfio.o modcdfnames.o $(FFLAGS)

cdfcurl: cdfio.o  cdfcurl.f90
	$(F90) cdfcurl.f90 -o $(BINDIR)/cdfcurl cdfio.o modcdfnames.o $(FFLAGS)

cdfdiv: cdfio.o  cdfdiv.f90
	$(F90) cdfdiv.f90 -o $(BINDIR)/cdfdiv cdfio.o modcdfnames.o $(FFLAGS)

cdflap: cdfio.o  cdflap.f90
	$(F90) cdflap.f90 -o $(BINDIR)/cdflap cdfio.o modcdfnames.o $(FFLAGS)

cdfeddyscale_pass1: cdfio.o  cdfeddyscale_pass1.f90
	$(F90) cdfeddyscale_pass1.f90 -o $(BINDIR)/cdfeddyscale_pass1 cdfio.o modcdfnames.o $(FFLAGS)

cdfeddyscale: cdfio.o  cdfeddyscale.f90
	$(F90) cdfeddyscale.f90 -o $(BINDIR)/cdfeddyscale cdfio.o modcdfnames.o $(FFLAGS)

cdfw: cdfio.o  cdfw.f90
	$(F90) cdfw.f90 -o $(BINDIR)/cdfw cdfio.o modcdfnames.o $(FFLAGS)

cdfgeo-uv: cdfio.o  cdfgeo-uv.f90
	$(F90) cdfgeo-uv.f90 -o $(BINDIR)/cdfgeo-uv cdfio.o modcdfnames.o $(FFLAGS)

cdfgeostrophy: cdfio.o eos.o cdfgeostrophy.f90
	$(F90) $(OMP) cdfgeostrophy.f90 -o $(BINDIR)/cdfgeostrophy cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfmxl: cdfio.o eos.o  cdfmxl.f90
	$(F90) $(OMP) cdfmxl.f90 -o $(BINDIR)/cdfmxl cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfrhoproj: cdfio.o  cdfrhoproj.f90
	$(F90) $(OMP) cdfrhoproj.f90 -o $(BINDIR)/cdfrhoproj cdfio.o modcdfnames.o  $(FFLAGS) 

cdfzisot: cdfio.o  cdfzisot.f90
	$(F90) cdfzisot.f90 -o $(BINDIR)/cdfzisot cdfio.o modcdfnames.o  $(FFLAGS) 

cdfsigintegr: cdfio.o  modutils.o cdfsigintegr.f90
	$(F90) cdfsigintegr.f90 -o $(BINDIR)/cdfsigintegr cdfio.o modcdfnames.o modutils.o  $(FFLAGS) 

cdfisopsi: cdfio.o eos.o cdfisopsi.f90
	$(F90) $(OMP) cdfisopsi.f90 -o $(BINDIR)/cdfisopsi cdfio.o eos.o modcdfnames.o  $(FFLAGS) 

cdfpvor: eos.o cdfio.o  cdfpvor.f90
	$(F90) $(OMP) cdfpvor.f90 -o $(BINDIR)/cdfpvor cdfio.o eos.o modcdfnames.o  $(FFLAGS) 

cdfgradT: cdfio.o  cdfgradT.f90
	$(F90) cdfgradT.f90 -o $(BINDIR)/cdfgradT cdfio.o  modcdfnames.o  $(FFLAGS) 

cdfkempemekeepe: cdfio.o  cdfkempemekeepe.f90
	$(F90) cdfkempemekeepe.f90 -o $(BINDIR)/cdfkempemekeepe cdfio.o  modcdfnames.o $(FFLAGS) 

cdfbci: cdfio.o  cdfbci.f90
	$(F90) cdfbci.f90 -o $(BINDIR)/cdfbci cdfio.o modcdfnames.o  $(FFLAGS) 

cdfbti: cdfio.o  cdfbti.f90
	$(F90) cdfbti.f90 -o $(BINDIR)/cdfbti cdfio.o modcdfnames.o $(FFLAGS) 

cdfnrjcomp: cdfio.o  cdfnrjcomp.f90
	$(F90) cdfnrjcomp.f90 -o $(BINDIR)/cdfnrjcomp cdfio.o modcdfnames.o $(FFLAGS) 

cdfhdy: cdfio.o  eos.o  cdfhdy.f90
	$(F90) $(OMP) cdfhdy.f90 -o $(BINDIR)/cdfhdy cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfhdy3d: cdfio.o  eos.o  cdfhdy3d.f90
	$(F90) $(OMP) cdfhdy3d.f90 -o $(BINDIR)/cdfhdy3d cdfio.o eos.o  modcdfnames.o $(FFLAGS)

cdfmaskdmp: cdfio.o  eos.o  cdfmaskdmp.f90
	$(F90) $(OMP) cdfmaskdmp.f90 -o $(BINDIR)/cdfmaskdmp cdfio.o eos.o modcdfnames.o  $(FFLAGS)

## Transport programs
cdfmhst: cdfio.o  cdfmhst.f90
	$(F90) cdfmhst.f90 -o $(BINDIR)/cdfmhst cdfio.o modcdfnames.o $(FFLAGS)

cdfvhst: cdfio.o  cdfvhst.f90
	$(F90) cdfvhst.f90 -o $(BINDIR)/cdfvhst cdfio.o modcdfnames.o $(FFLAGS)

cdfvtrp: cdfio.o  cdfvtrp.f90
	$(F90) cdfvtrp.f90 -o $(BINDIR)/cdfvtrp cdfio.o modcdfnames.o $(FFLAGS)

cdfpsi: cdfio.o  modutils.o cdfpsi.f90
	$(F90) cdfpsi.f90  -o $(BINDIR)/cdfpsi cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfpsi_level: cdfio.o  cdfpsi_level.f90
	$(F90) cdfpsi_level.f90  -o $(BINDIR)/cdfpsi_level cdfio.o modcdfnames.o $(FFLAGS)

cdftransport: cdfio.o  modutils.o cdftransport.f90
	$(F90) cdftransport.f90 -o $(BINDIR)/cdftransport cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfvFWov: cdfio.o  modutils.o cdfvFWov.f90
	$(F90) cdfvFWov.f90 -o $(BINDIR)/cdfvFWov cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfsigtrp: cdfio.o eos.o  modutils.o cdfsigtrp.f90 
	$(F90)  $(OMP) cdfsigtrp.f90 -o $(BINDIR)/cdfsigtrp cdfio.o eos.o modcdfnames.o modutils.o  $(FFLAGS)

cdfsigtrp_broken: cdfio.o eos.o  modutils.o cdfsigtrp_broken.f90 
	$(F90)  $(OMP) cdfsigtrp_broken.f90 -o $(BINDIR)/cdfsigtrp_broken cdfio.o eos.o modcdfnames.o modutils.o  $(FFLAGS)

cdftransig_xy3d: cdfio.o eos.o  modutils.o cdftransig_xy3d.f90
	$(F90)  $(OMP) cdftransig_xy3d.f90 -o $(BINDIR)/cdftransig_xy3d cdfio.o eos.o modcdfnames.o modutils.o  $(FFLAGS)

cdftempvol-full: cdfio.o  cdftempvol-full.f90
	$(F90) cdftempvol-full.f90 -o $(BINDIR)/cdftempvol-full cdfio.o modcdfnames.o  $(FFLAGS)

cdfmoc: cdfio.o  eos.o cdftools.o cdfmoc.f90
	$(F90) $(OMP) cdfmoc.f90 -o $(BINDIR)/cdfmoc cdfio.o eos.o modcdfnames.o cdftools.o $(FFLAGS)

cdfmht_gsop: cdfio.o eos.o  cdfmht_gsop.f90
	$(F90) $(OMP) cdfmht_gsop.f90 -o $(BINDIR)/cdfmht_gsop cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfmoc_rapid_26N_r8_ORCA025: cdfio.o eos.o  cdfmoc_rapid_26N_r8_ORCA025.f90
	$(F90) $(OMP) cdfmoc_rapid_26N_r8_ORCA025.f90 -o $(BINDIR)/cdfmoc_rapid_26N_r8_ORCA025 cdfio.o eos.o $(FFLAGS)

cdfmocsig: cdfio.o eos.o modutils.o  cdfmocsig.f90 
	$(F90) $(OMP) cdfmocsig.f90 -o $(BINDIR)/cdfmocsig cdfio.o eos.o modcdfnames.o modutils.o $(FFLAGS)

cdfmean: cdfio.o  cdfmean.f90
	$(F90) cdfmean.f90 -o $(BINDIR)/cdfmean cdfio.o modcdfnames.o $(FFLAGS)

cdfdegradt: cdfio.o  cdfdegradt.f90
	$(F90) cdfdegradt.f90 -o $(BINDIR)/cdfdegradt cdfio.o modcdfnames.o $(FFLAGS)

cdfdegradu: cdfio.o  cdfdegradu.f90
	$(F90) cdfdegradu.f90 -o $(BINDIR)/cdfdegradu cdfio.o modcdfnames.o $(FFLAGS)

cdfdegradv: cdfio.o  cdfdegradv.f90
	$(F90) cdfdegradv.f90 -o $(BINDIR)/cdfdegradv cdfio.o modcdfnames.o $(FFLAGS)

cdfdegradw: cdfio.o  cdfdegradw.f90
	$(F90) cdfdegradw.f90 -o $(BINDIR)/cdfdegradw cdfio.o modcdfnames.o $(FFLAGS)

cdfsum: cdfio.o  cdfsum.f90
	$(F90) cdfsum.f90 -o $(BINDIR)/cdfsum cdfio.o modcdfnames.o $(FFLAGS)

cdfvertmean: cdfio.o  modutils.o cdfvertmean.f90
	$(F90) cdfvertmean.f90 -o $(BINDIR)/cdfvertmean cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfvint: cdfio.o  modutils.o cdfvint.f90
	$(F90) cdfvint.f90 -o $(BINDIR)/cdfvint cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfheatc: cdfio.o modutils.o cdfheatc.f90
	$(F90) cdfheatc.f90 -o $(BINDIR)/cdfheatc cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfmxlheatc: cdfio.o  modutils.o cdfmxlheatc.f90
	$(F90) cdfmxlheatc.f90 -o $(BINDIR)/cdfmxlheatc cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfmxlsaltc: cdfio.o  modutils.o cdfmxlsaltc.f90
	$(F90) cdfmxlsaltc.f90 -o $(BINDIR)/cdfmxlsaltc cdfio.o modcdfnames.o modutils.o $(FFLAGS)

cdfmxlhcsc: cdfio.o  eos.o cdfmxlhcsc.f90
	$(F90) $(OMP) cdfmxlhcsc.f90 -o $(BINDIR)/cdfmxlhcsc cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdficediags: cdfio.o  cdficediags.f90
	$(F90) cdficediags.f90 -o $(BINDIR)/cdficediags cdfio.o modcdfnames.o $(FFLAGS)

cdfzonalmean: cdfio.o  cdfzonalmean.f90
	$(F90) $(OMP) cdfzonalmean.f90 -o $(BINDIR)/cdfzonalmean cdfio.o modcdfnames.o $(FFLAGS) 

cdfzonalmeanvT: cdfio.o  modutils.o  cdfzonalmeanvT.f90
	$(F90) cdfzonalmeanvT.f90 -o $(BINDIR)/cdfzonalmeanvT cdfio.o modcdfnames.o modutils.o $(FFLAGS) 

cdfzonalmeanvT2: cdfio.o  modutils.o  cdfzonalmeanvT2.f90
	$(F90) cdfzonalmeanvT2.f90 -o $(BINDIR)/cdfzonalmeanvT2 cdfio.o modcdfnames.o modutils.o $(FFLAGS) 

cdfzonalsum: cdfio.o  cdfzonalsum.f90
	$(F90) cdfzonalsum.f90 -o $(BINDIR)/cdfzonalsum cdfio.o modcdfnames.o $(FFLAGS) 

cdfzonalout: cdfio.o  cdfzonalout.f90
	$(F90) cdfzonalout.f90 -o $(BINDIR)/cdfzonalout cdfio.o modcdfnames.o $(FFLAGS) 

cdfhflx: cdfio.o  cdfhflx.f90
	$(F90) cdfhflx.f90 -o $(BINDIR)/cdfhflx cdfio.o modcdfnames.o $(FFLAGS)

cdfwflx: cdfio.o  cdfwflx.f90
	$(F90) cdfwflx.f90 -o $(BINDIR)/cdfwflx cdfio.o modcdfnames.o $(FFLAGS)

cdfbuoyflx: cdfio.o  eos.o cdfbuoyflx.f90
	$(F90) $(OMP) cdfbuoyflx.f90 -o $(BINDIR)/cdfbuoyflx cdfio.o eos.o modcdfnames.o $(FFLAGS)

## Extracting tools, information tools
cdfprofile: cdfio.o  cdfprofile.f90
	$(F90) cdfprofile.f90  -o $(BINDIR)/cdfprofile cdfio.o modcdfnames.o $(FFLAGS)

cdfwhereij:cdfio.o  cdfwhereij.f90
	$(F90) cdfwhereij.f90  -o $(BINDIR)/cdfwhereij cdfio.o modcdfnames.o $(FFLAGS)

cdffindij: cdfio.o cdftools.o  cdffindij.f90
	$(F90) cdffindij.f90  -o $(BINDIR)/cdffindij cdfio.o cdftools.o modcdfnames.o $(FFLAGS)

cdf_use_lib: cdftools.o cdf_use_lib.f90
	$(F90) cdf_use_lib.f90  -o $(BINDIR)/cdf_use_lib cdfio.o  cdftools.o $(FFLAGS)

cdfweight: cdfio.o  cdftools.o cdfweight.f90 
	$(F90) cdfweight.f90  -o $(BINDIR)/cdfweight cdfio.o cdftools.o modcdfnames.o $(FFLAGS)

cdfweight2D: cdfio.o  cdfweight2D.f90
	$(F90) cdfweight2D.f90  -o $(BINDIR)/cdfweight2D cdfio.o $(FFLAGS)

cdfcoloc: cdfio.o  cdfcoloc.f90
	$(F90) cdfcoloc.f90  -o $(BINDIR)/cdfcoloc cdfio.o modcdfnames.o $(FFLAGS)

cdfcoloc2D: cdfio.o  cdfcoloc2D.f90
	$(F90)  cdfcoloc2D.f90  -o $(BINDIR)/cdfcoloc2D cdfio.o $(FFLAGS)

cdfcoloc2: cdfio.o  cdfcoloc2.f90
	$(F90) cdfcoloc2.f90  -o $(BINDIR)/cdfcoloc2 cdfio.o $(FFLAGS)

cdfcoloc3: cdfio.o  cdfcoloc3.f90
	$(F90) cdfcoloc3.f90  -o $(BINDIR)/cdfcoloc3 cdfio.o $(FFLAGS)

cdf2levitusgrid2d: cdfio.o cdftools.o modutils.o cdf2levitusgrid2d.f90
	$(F90) cdf2levitusgrid2d.f90  -o $(BINDIR)/cdf2levitusgrid2d  cdfio.o modcdfnames.o cdftools.o modutils.o $(FFLAGS)

cdf2levitusgrid3d: cdfio.o cdftools.o modutils.o cdf2levitusgrid3d.f90
	$(F90) cdf2levitusgrid3d.f90  -o $(BINDIR)/cdf2levitusgrid3d  cdfio.o modcdfnames.o cdftools.o modutils.o $(FFLAGS)

cdf025to05: cdfio.o  cdf025to05.f90
	$(F90) cdf025to05.f90  -o $(BINDIR)/cdf025to05 cdfio.o modcdfnames.o $(FFLAGS)

cdfrnf025to05: cdfio.o  cdfrnf025to05.f90
	$(F90) cdfrnf025to05.f90  -o $(BINDIR)/cdfrnf025to05 cdfio.o modcdfnames.o $(FFLAGS)

cdfstatcoord: cdfio.o  cdfstatcoord.f90
	$(F90) cdfstatcoord.f90  -o $(BINDIR)/cdfstatcoord cdfio.o modcdfnames.o $(FFLAGS)

cdfmaxmoc: cdfio.o cdfmaxmoc.f90 
	$(F90) cdfmaxmoc.f90  -o $(BINDIR)/cdfmaxmoc cdfio.o modcdfnames.o $(FFLAGS)

cdfcensus: cdfio.o eos.o cdfcensus.f90
	$(F90) $(OMP) cdfcensus.f90  -o $(BINDIR)/cdfcensus cdfio.o eos.o modcdfnames.o $(FFLAGS)

cdfzoom: cdfio.o  cdfzoom.f90
	$(F90) cdfzoom.f90  -o $(BINDIR)/cdfzoom cdfio.o modcdfnames.o $(FFLAGS)

cdfmax: cdfio.o  cdfmax.f90
	$(F90) cdfmax.f90  -o $(BINDIR)/cdfmax cdfio.o modcdfnames.o $(FFLAGS)

cdfprobe: cdfio.o  cdfprobe.f90
	$(F90) cdfprobe.f90  -o $(BINDIR)/cdfprobe cdfio.o modcdfnames.o $(FFLAGS)

cdfinfo: cdfio.o  cdfinfo.f90
	$(F90) cdfinfo.f90  -o $(BINDIR)/cdfinfo cdfio.o modcdfnames.o $(FFLAGS)

cdfclip: cdfio.o  cdfclip.f90
	$(F90) cdfclip.f90  -o $(BINDIR)/cdfclip cdfio.o modcdfnames.o $(FFLAGS)

cdfsmooth: cdfio.o  modutils.o cdfsmooth.f90
	$(F90) cdfsmooth.f90  -o $(BINDIR)/cdfsmooth cdfio.o modutils.o modcdfnames.o $(FFLAGS)

cdfpendep: cdfio.o  cdfpendep.f90
	$(F90) cdfpendep.f90  -o $(BINDIR)/cdfpendep cdfio.o modcdfnames.o $(FFLAGS)

cdffracinv: cdfio.o cdffracinv.f90
	$(F90) cdffracinv.f90 -o $(BINDIR)/cdffracinv cdfio.o modcdfnames.o $(FFLAGS)

cdfzgrv3: cdfio.o  cdfzgrv3.f90
	$(F90) cdfzgrv3.f90  -o $(BINDIR)/cdfzgrv3 cdfio.o $(FFLAGS)

## reformating programs
cdf16bit: cdfio.o cdf16bit.f90
	$(F90) cdf16bit.f90  -o $(BINDIR)/cdf16bit cdfio.o modcdfnames.o $(FFLAGS)

cdf2matlab: cdfio.o cdf2matlab.f90
	$(F90)    cdf2matlab.f90  -o $(BINDIR)/cdf2matlab cdfio.o modcdfnames.o $(FFLAGS)

cdfvita: cdfio.o cdfvita.f90
	$(F90) cdfvita.f90  -o $(BINDIR)/cdfvita cdfio.o modcdfnames.o $(FFLAGS)

cdfvita-geo: cdfio.o cdfvita-geo.f90
	$(F90) cdfvita-geo.f90  -o $(BINDIR)/cdfvita-geo  cdfio.o modcdfnames.o $(FFLAGS)

cdfconvert: cdfio.o cdfconvert.f90
	$(F90)  cdfconvert.f90  -o $(BINDIR)/cdfconvert cdfio.o modcdfnames.o $(FFLAGS)

cdfflxconv: cdfio.o cdfflxconv.f90
	$(F90)   cdfflxconv.f90  -o $(BINDIR)/cdfflxconv cdfio.o modcdfnames.o $(FFLAGS)

cdfsstconv: cdfio.o cdfsstconv.f90
	$(F90)   cdfsstconv.f90  -o $(BINDIR)/cdfsstconv cdfio.o modcdfnames.o $(FFLAGS)

cdfstrconv: cdfio.o cdfstrconv.f90
	$(F90)   cdfstrconv.f90  -o $(BINDIR)/cdfstrconv cdfio.o modcdfnames.o $(FFLAGS)

cdfbathy: cdfio.o cdfbathy.f90
	$(F90)   cdfbathy.f90  -o $(BINDIR)/cdfbathy cdfio.o modcdfnames.o $(FFLAGS)

cdfcofdis: cdfio.o cdfcofdis.f90
	$(F90)    cdfcofdis.f90  -o $(BINDIR)/cdfcofdis cdfio.o modcdfnames.o $(FFLAGS)

cdfcoastline: cdfio.o cdfcoastline.f90
	$(F90)    cdfcoastline.f90  -o $(BINDIR)/cdfcoastline cdfio.o modcdfnames.o $(FFLAGS)

cdfvar: cdfbathy
	ln -sf  cdfbathy $(BINDIR)/cdfvar 

cdfcsp: cdfio.o cdfcsp.f90
	$(F90)   cdfcsp.f90  -o $(BINDIR)/cdfcsp cdfio.o modcdfnames.o $(FFLAGS)

cdfnan: cdfio.o cdfnan.f90
	$(F90)   cdfnan.f90  -o $(BINDIR)/cdfnan cdfio.o modcdfnames.o $(FFLAGS)

cdfscale: cdfio.o cdfscale.f90
	$(F90)   cdfscale.f90  -o $(BINDIR)/cdfscale cdfio.o modcdfnames.o $(FFLAGS)

cdfnorth_unfold: cdfio.o cdfnorth_unfold.f90
	$(F90)   cdfnorth_unfold.f90  -o $(BINDIR)/cdfnorth_unfold cdfio.o modcdfnames.o $(FFLAGS)

cdfpolymask: cdfio.o modpoly.o cdfpolymask.f90
	$(F90)   cdfpolymask.f90  -o $(BINDIR)/cdfpolymask cdfio.o modpoly.o modcdfnames.o $(FFLAGS)

cdfovide: cdfio.o  cdfovide.f90
	$(F90) cdfovide.f90  -o $(BINDIR)/cdfovide cdfio.o modcdfnames.o $(FFLAGS)

cdf_xtrac_brokenline: cdfio.o  modcdfnames.o cdftools.o  cdf_xtrac_brokenline.f90
	$(F90)  cdf_xtrac_brokenline.f90  -o $(BINDIR)/cdf_xtrac_brokenline cdfio.o cdftools.o modcdfnames.o $(FFLAGS)

cdfmppini: cdfio.o  cdfmppini.f90
	$(F90)  cdfmppini.f90  -o $(BINDIR)/cdfmppini cdfio.o modcdfnames.o $(FFLAGS)

cdffixtime: cdfio.o  cdffixtime.f90
	$(F90)  cdffixtime.f90  -o $(BINDIR)/cdffixtime cdfio.o modcdfnames.o $(FFLAGS)

cdfnamelist: modcdfnames.o  cdfnamelist.f90
	$(F90)  cdfnamelist.f90  -o $(BINDIR)/cdfnamelist  modcdfnames.o $(FFLAGS) $(FDATE_FLAG)

cdfchgrid: cdfio.o cdfchgrid.f90
	$(F90) cdfchgrid.f90  -o $(BINDIR)/cdfchgrid cdfio.o modcdfnames.o $(FFLAGS)

cdfcmp: cdfio.o cdfcmp.f90
	$(F90) cdfcmp.f90  -o $(BINDIR)/cdfcmp cdfio.o modcdfnames.o $(FFLAGS)

cdfsections: eos.o cdfsections.f90
	$(F90) $(OMP) cdfsections.f90  -o $(BINDIR)/cdfsections eos.o modcdfnames.o $(FFLAGS)

## Modules

cdfio.o: cdfio.F90 modcdfnames.o
	$(F90) -c  cdfio.F90 $(FFLAGS)

eos.o: eos.f90
	$(F90) $(OMP) -c eos.f90 $(FFLAGS)

cdftools.o: cdfio.o cdftools.f90
	$(F90) -c cdftools.f90 $(FFLAGS)

modpoly.o: modpoly.f90
	$(F90) -c modpoly.f90 $(FFLAGS)

modcdfnames.o: modcdfnames.f90
	$(F90) -c modcdfnames.f90 $(FFLAGS)

modutils.o: cdfio.o modutils.f90
	$(F90) -c modutils.f90 $(FFLAGS)

## Utilities
clean:
	\rm -f *.mod *.o  *~ *.1 *.opod

cleanexe: clean
	( cd $(BINDIR) ; \rm -f $(EXEC) )

man: cdftools.1

cdftools.1: cdftools.opod
	pod2man --center "CDFTOOLS / NEMO Documentation" \
	  --release "SVN Revision $$(LANG=C svn update | grep '^At rev' | awk '{print $$3}' | cut -f 1 -d '.')" \
	  cdftools.opod > cdftools.1

cdftools.opod: $(EXEC) cdftools-begin.pod cdftools-end.pod
	cat cdftools-begin.pod > cdftools.opod
	for s in $$( cd $(BINDIR); ls -1 ); do echo ''; echo "=head2 $$s"; echo ''; $$s; done >> cdftools.opod
	cat cdftools-end.pod >> cdftools.opod

install:
	@mkdir -p $(INSTALL)
	cd bin ; \cp $(EXEC)  $(INSTALL)

installman:
	@mkdir -p $(INSTALL_MAN)/man1;
	\cp -f cdftools.1 $(INSTALL_MAN)/man1/;
	for s in $$( cd $(BINDIR); ls -1 ); do ( cd $(INSTALL_MAN)/man1/; ln -sf cdftools.1 $$s.1 ); done;
