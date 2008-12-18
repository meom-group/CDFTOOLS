# Makefile for CDFTOOLS 

# ( make.macro is a link that points to the file macro.xxx where 
#   xxx is representative of your machine )
# !!  $Rev$
# !!  $Date$
# !!  $Id$
# !!--------------------------------------------------------------


include make.macro


CDFTOOLS=CDFTOOLS-2.1

EXEC = cdfmoy cdfmoyt  cdfmoy_sp cdfstd cdfmoy_sal2_temp2  cdfmoy_annual cdfmoy_chsp cdfmoy_freq cdfvT cdfvsig cdfspeed cdfsum\
       cdfmoyuv cdfmoyuvwt \
       cdfeke cdfrmsssh cdfstdevw cdfstdevts cdflinreg cdfimprovechk\
       cdfbn2 cdfbn2-full  cdfsig0 cdfsigi cdfsiginsitu cdfbottomsig0 cdfbottomsigi cdfbottom cdfets cdfcurl cdfw cdfgeo-uv cdfmxl cdfmxl-full\
       cdfrhoproj cdfisopycdep cdfsigintegr cdfpv cdflspv cdfpvor\
       cdfmhst cdfmhst-full cdfvhst cdfvhst-full cdftransportiz cdftransportiz_noheat cdftransportiz-full \
       cdftransportizpm cdfmoc_gsop  cdfmoc_gsop_x cdfmht_gsop \
       cdfmasstrp cdfmasstrp-full \
       cdfsigtrp cdfsigitrp cdfsigtrp-full cdftemptrp-full  cdftempvol-full\
       cdfpsi cdfpsi-full cdfpsi-open cdfmoc cdfmoc-full cdfmocatl cdfmocsig cdfmean cdfmeanvar cdfmean-full cdfzeromean \
       cdfheatc cdfheatc-full cdfzonalmean cdfhflx cdfwflx cdfbuoyflx\
       cdfmxlheatc cdfmxlheatc-full cdfmxlsaltc cdfmxlhcsc cdfvertmean\
       cdfpendep cdfzonalsum cdficediags cdfzonalout\
       cdfprofile  cdfwhereij cdffindij cdfweight cdfmaxmoc cdfcensus cdfzoom cdfmax cdfmax_sp cdfprobe \
       bimgmoy4 bimgcaltrans cdf16bit cdfvita cdfconvert cdfflxconv cdfclip cdfsstconv cdfstrconv cdfbathy cdfvar cdfmkmask-zone\
       cdfcsp cdfcoloc cdfmltmask cdfstatcoord  cdfpolymask cdfsmooth cdfmkmask \
       cdfkempemekeepe cdfbci cdfbti cdfnrjcomp cdfcofdis

all: $(EXEC)

## Statistical programs
cdfmoy: cdfio.o   cdfmoy.f90
	$(F90) cdfmoy.f90 -o cdfmoy cdfio.o  $(FFLAGS)

cdfmoyt: cdfio.o   cdfmoyt.f90
	$(F90) cdfmoyt.f90 -o cdfmoyt cdfio.o  $(FFLAGS)

cdfmoy_mpp: cdfio.o   cdfmoy_mpp.f90
	$(MPF90) cdfmoy_mpp.f90 -o cdfmoy_mpp cdfio.o  $(FFLAGS) $(LMPI)

cdfmoy_sal2_temp2: cdfio.o   cdfmoy_sal2_temp2.f90
	$(F90) cdfmoy_sal2_temp2.f90 -o cdfmoy_sal2_temp2 cdfio.o  $(FFLAGS)

cdfmoy_sp: cdfio.o   cdfmoy_sp.f90
	$(F90) cdfmoy_sp.f90 -o cdfmoy_sp cdfio.o  $(FFLAGS)

cdfmoy_chsp: cdfio.o   cdfmoy_chsp.f90
	$(F90) cdfmoy_chsp.f90 -o cdfmoy_chsp cdfio.o  $(FFLAGS)

cdfmoy_freq: cdfio.o   cdfmoy_freq.f90
	$(F90) cdfmoy_freq.f90 -o cdfmoy_freq cdfio.o  $(FFLAGS)

cdfmoyuv: cdfio.o   cdfmoyuv.f90
	$(F90) cdfmoyuv.f90 -o cdfmoyuv cdfio.o  $(FFLAGS)

cdfmoyuvwt: cdfio.o   cdfmoyuvwt.f90
	$(F90) cdfmoyuvwt.f90 -o cdfmoyuvwt cdfio.o  $(FFLAGS)

cdfstd: cdfio.o  cdfstd.f90
	$(F90)  cdfstd.f90 -o cdfstd cdfio.o $(FFLAGS)

cdfmoy_annual: cdfio.o   cdfmoy_annual.f90
	$(F90) cdfmoy_annual.f90 -o cdfmoy_annual cdfio.o  $(FFLAGS)

cdfeke: cdfio.o  cdfeke.f90
	$(F90) cdfeke.f90 -o cdfeke cdfio.o $(FFLAGS)

cdfrmsssh: cdfio.o  cdfrmsssh.f90
	$(F90) cdfrmsssh.f90 -o cdfrmsssh cdfio.o $(FFLAGS)

cdfstdevw: cdfio.o  cdfstdevw.f90
	$(F90) cdfstdevw.f90 -o cdfstdevw cdfio.o $(FFLAGS)

cdfstdevts: cdfio.o  cdfstdevts.f90
	$(F90) cdfstdevts.f90 -o cdfstdevts cdfio.o $(FFLAGS)

cdfvT: cdfio.o  cdfvT.f90
	$(F90) cdfvT.f90 -o cdfvT cdfio.o $(FFLAGS)

cdfvsig: cdfio.o eos.o  cdfvsig.f90
	$(F90) cdfvsig.f90 -o cdfvsig cdfio.o eos.o $(FFLAGS)

cdfspeed: cdfio.o  cdfspeed.f90
	$(F90) cdfspeed.f90 -o cdfspeed cdfio.o $(FFLAGS)

cdfimprovechk: cdfio.o  cdfimprovechk.f90
	$(F90) cdfimprovechk.f90 -o cdfimprovechk cdfio.o $(FFLAGS)

cdflinreg: cdfio.o  cdflinreg.f90
	$(F90) cdflinreg.f90 -o cdflinreg cdfio.o $(FFLAGS)

## Derived quantities programs
cdfbn2: cdfio.o  eos.o  cdfbn2.f90
	$(F90) cdfbn2.f90 -o cdfbn2 cdfio.o eos.o  $(FFLAGS)

cdfbn2-full: cdfio.o  eos.o  cdfbn2-full.f90
	$(F90) cdfbn2-full.f90 -o cdfbn2-full cdfio.o eos.o  $(FFLAGS)

cdfsig0: cdfio.o  eos.o  cdfsig0.f90
	$(F90) cdfsig0.f90 -o cdfsig0 cdfio.o eos.o  $(FFLAGS)

cdfsigi: cdfio.o  eos.o  cdfsigi.f90
	$(F90) cdfsigi.f90 -o cdfsigi cdfio.o eos.o  $(FFLAGS)

cdfsiginsitu: cdfio.o  eos.o  cdfsiginsitu.f90
	$(F90) cdfsiginsitu.f90 -o cdfsiginsitu cdfio.o eos.o  $(FFLAGS)

cdfbottomsig0: cdfio.o  eos.o  cdfbottomsig0.f90
	$(F90) cdfbottomsig0.f90 -o cdfbottomsig0 cdfio.o eos.o  $(FFLAGS)

cdfbottomsigi: cdfio.o  eos.o  cdfbottomsigi.f90
	$(F90) cdfbottomsigi.f90 -o cdfbottomsigi cdfio.o eos.o  $(FFLAGS)

cdfbottom: cdfio.o    cdfbottom.f90
	$(F90) cdfbottom.f90 -o cdfbottom cdfio.o   $(FFLAGS)

cdfets: cdfio.o  eos.o  cdfets.f90
	$(F90) cdfets.f90 -o cdfets cdfio.o eos.o  $(FFLAGS)

cdfmsk: cdfio.o  cdfmsk.f90
	$(F90) cdfmsk.f90 -o cdfmsk cdfio.o $(FFLAGS)

cdfmsksal: cdfio.o  cdfmsksal.f90
	$(F90) cdfmsksal.f90 -o cdfmsksal cdfio.o $(FFLAGS)

cdfmkmask: cdfio.o  cdfmkmask.f90
	$(F90) cdfmkmask.f90 -o cdfmkmask cdfio.o $(FFLAGS)

cdfmkmask-zone: cdfio.o  cdfmkmask-zone.f90
	$(F90) cdfmkmask-zone.f90 -o cdfmkmask-zone cdfio.o $(FFLAGS)

cdfmltmask: cdfio.o  cdfmltmask.f90
	$(F90) cdfmltmask.f90 -o cdfmltmask cdfio.o $(FFLAGS)

cdfmltmask2: cdfio.o  cdfmltmask2.f90
	$(F90) cdfmltmask2.f90 -o cdfmltmask2 cdfio.o $(FFLAGS)

cdfcurl: cdfio.o  cdfcurl.f90
	$(F90) cdfcurl.f90 -o cdfcurl cdfio.o $(FFLAGS)

cdfw: cdfio.o  cdfw.f90
	$(F90) cdfw.f90 -o cdfw cdfio.o $(FFLAGS)

cdfgeo-uv: cdfio.o  cdfgeo-uv.f90
	$(F90) cdfgeo-uv.f90 -o cdfgeo-uv cdfio.o $(FFLAGS)

cdfmxl: cdfio.o eos.o  cdfmxl.f90
	$(F90) cdfmxl.f90 -o cdfmxl cdfio.o eos.o $(FFLAGS)

cdfmxl-full: cdfio.o eos.o  cdfmxl-full.f90
	$(F90) cdfmxl-full.f90 -o cdfmxl-full cdfio.o eos.o $(FFLAGS)

cdfrhoproj: cdfio.o  cdfrhoproj.f90
	$(F90) cdfrhoproj.f90 -o cdfrhoproj cdfio.o  $(FFLAGS) 

cdfisopycdep: cdfio.o  cdfisopycdep.f90
	$(F90) cdfisopycdep.f90 -o cdfisopycdep cdfio.o  $(FFLAGS) 

cdfsigintegr: cdfio.o  cdfsigintegr.f90
	$(F90) cdfsigintegr.f90 -o cdfsigintegr cdfio.o  $(FFLAGS) 

cdfpv: cdfio.o  cdfpv.f90
	$(F90) cdfpv.f90 -o cdfpv cdfio.o eos.o  $(FFLAGS) 

cdflspv: cdfio.o  cdflspv.f90
	$(F90) cdflspv.f90 -o cdflspv cdfio.o eos.o  $(FFLAGS) 

cdfpvor: cdfio.o  cdfpvor.f90
	$(F90) cdfpvor.f90 -o cdfpvor cdfio.o eos.o  $(FFLAGS) 

cdfkempemekeepe: cdfio.o  cdfkempemekeepe.f90
	$(F90) cdfkempemekeepe.f90 -o cdfkempemekeepe cdfio.o  $(FFLAGS) 

cdfbci: cdfio.o  cdfbci.f90
	$(F90) cdfbci.f90 -o cdfbci cdfio.o  $(FFLAGS) 

cdfbti: cdfio.o  cdfbti.f90
	$(F90) cdfbti.f90 -o cdfbti cdfio.o  $(FFLAGS) 

cdfnrjcomp: cdfio.o  cdfnrjcomp.f90
	$(F90) cdfnrjcomp.f90 -o cdfnrjcomp cdfio.o  $(FFLAGS) 

## Transport programs
cdfmhst: cdfio.o  cdfmhst.f90
	$(F90) cdfmhst.f90 -o cdfmhst cdfio.o $(FFLAGS)

cdfmhst-full: cdfio.o  cdfmhst-full.f90
	$(F90) cdfmhst-full.f90 -o cdfmhst-full cdfio.o $(FFLAGS)

cdfvhst: cdfio.o  cdfvhst.f90
	$(F90) cdfvhst.f90 -o cdfvhst cdfio.o $(FFLAGS)

cdfvtrp: cdfio.o  cdfvtrp.f90
	$(F90) cdfvtrp.f90 -o cdfvtrp cdfio.o $(FFLAGS)

cdfvhst-full: cdfio.o  cdfvhst-full.f90
	$(F90) cdfvhst-full.f90 -o cdfvhst-full cdfio.o $(FFLAGS)

cdfpsi: cdfio.o  cdfpsi.f90
	$(F90) cdfpsi.f90  -o cdfpsi cdfio.o $(FFLAGS)

cdfpsi-full: cdfio.o  cdfpsi-full.f90
	$(F90) cdfpsi-full.f90  -o cdfpsi-full cdfio.o $(FFLAGS)

cdfpsi-open: cdfio.o  cdfpsi-open.f90
	$(F90) cdfpsi-open.f90  -o cdfpsi-open cdfio.o $(FFLAGS)

cdfpsi-open-zap: cdfio.o  cdfpsi-open-zap.f90
	$(F90) cdfpsi-open-zap.f90  -o cdfpsi-open-zap cdfio.o $(FFLAGS)

cdfpsi-open_AM: cdfio.o  cdfpsi-open_AM.f90
	$(F90) cdfpsi-open_AM.f90  -o cdfpsi-open_AM cdfio.o $(FFLAGS)

cdftransportiz: cdfio.o  cdftransportiz.f90
	$(F90) cdftransportiz.f90 -o cdftransportiz cdfio.o $(FFLAGS)

cdftransportizpm: cdfio.o  cdftransportizpm.f90
	$(F90) cdftransportizpm.f90 -o cdftransportizpm cdfio.o $(FFLAGS)

cdftransportiz_noheat: cdfio.o  cdftransportiz_noheat.f90
	$(F90) cdftransportiz_noheat.f90 -o cdftransportiz_noheat cdfio.o $(FFLAGS)

cdftransportiz-full: cdfio.o  cdftransportiz-full.f90
	$(F90) cdftransportiz-full.f90 -o cdftransportiz-full cdfio.o $(FFLAGS)

cdfmasstrp: cdfio.o  cdfmasstrp.f90
	$(F90) cdfmasstrp.f90 -o cdfmasstrp cdfio.o $(FFLAGS)

cdfmasstrp-full: cdfio.o  cdfmasstrp-full.f90
	$(F90) cdfmasstrp-full.f90 -o cdfmasstrp-full cdfio.o $(FFLAGS)

cdfmasstrp-julien: cdfio.o  cdfmasstrp-julien.f90
	$(F90) cdfmasstrp-julien.f90 -o cdfmasstrp-julien cdfio.o $(FFLAGS)

cdfsigtrp: cdfio.o  cdfsigtrp.f90
	$(F90) cdfsigtrp.f90 -o cdfsigtrp cdfio.o eos.o $(FFLAGS)

cdfsigitrp: cdfio.o  cdfsigitrp.f90
	$(F90) cdfsigitrp.f90 -o cdfsigitrp cdfio.o eos.o $(FFLAGS)

cdfsigtrp-full: cdfio.o  cdfsigtrp-full.f90
	$(F90) cdfsigtrp-full.f90 -o cdfsigtrp-full cdfio.o eos.o $(FFLAGS)

cdftemptrp-full: cdfio.o  cdftemptrp-full.f90
	$(F90) cdftemptrp-full.f90 -o cdftemptrp-full cdfio.o  $(FFLAGS)

cdftempvol-full: cdfio.o  cdftempvol-full.f90
	$(F90) cdftempvol-full.f90 -o cdftempvol-full cdfio.o  $(FFLAGS)

cdfmoc: cdfio.o  cdfmoc.f90
	$(F90) cdfmoc.f90 -o cdfmoc cdfio.o $(FFLAGS)

cdfmoc_gsop: cdfio.o eos.o  cdfmoc_gsop.f90
	$(F90) cdfmoc_gsop.f90 -o cdfmoc_gsop cdfio.o eos.o $(FFLAGS)

cdfmht_gsop: cdfio.o eos.o  cdfmht_gsop.f90
	$(F90) cdfmht_gsop.f90 -o cdfmht_gsop cdfio.o eos.o $(FFLAGS)

cdfmoc_gsop_x: cdfio.o eos.o  cdfmoc_gsop_x.f90
	$(F90) cdfmoc_gsop_x.f90 -o cdfmoc_gsop_x cdfio.o eos.o $(FFLAGS)

cdfmoc_rapid_26N_r8_ORCA025: cdfio.o eos.o  cdfmoc_rapid_26N_r8_ORCA025.f90
	$(F90) cdfmoc_rapid_26N_r8_ORCA025.f90 -o cdfmoc_rapid_26N_r8_ORCA025 cdfio.o eos.o $(FFLAGS)

cdfmocsig: cdfio.o eos.o  cdfmocsig.f90
	$(F90) cdfmocsig.f90 -o cdfmocsig cdfio.o eos.o $(FFLAGS)

cdfmoc-full: cdfio.o  cdfmoc-full.f90
	$(F90) cdfmoc-full.f90 -o cdfmoc-full cdfio.o $(FFLAGS)

cdfmocatl: cdfio.o  cdfmocatl.f90
	$(F90) cdfmocatl.f90 -o cdfmocatl cdfio.o $(FFLAGS)

cdfmean: cdfio.o  cdfmean.f90
	$(F90) cdfmean.f90 -o cdfmean cdfio.o $(FFLAGS)

cdfsum: cdfio.o  cdfsum.f90
	$(F90) cdfsum.f90 -o cdfsum cdfio.o $(FFLAGS)

cdfvertmean: cdfio.o  cdfvertmean.f90
	$(F90) cdfvertmean.f90 -o cdfvertmean cdfio.o $(FFLAGS)

cdfmeanvar: cdfio.o  cdfmeanvar.f90
	$(F90) cdfmeanvar.f90 -o cdfmeanvar cdfio.o $(FFLAGS)

cdfmean-full: cdfio.o  cdfmean-full.f90
	$(F90) cdfmean-full.f90 -o cdfmean-full cdfio.o $(FFLAGS)

cdfzeromean: cdfio.o  cdfzeromean.f90
	$(F90) cdfzeromean.f90 -o cdfzeromean cdfio.o $(FFLAGS)

cdfheatc: cdfio.o  cdfheatc.f90
	$(F90) cdfheatc.f90 -o cdfheatc cdfio.o $(FFLAGS)

cdfheatc-full: cdfio.o  cdfheatc-full.f90
	$(F90) cdfheatc-full.f90 -o cdfheatc-full cdfio.o $(FFLAGS)

cdfmxlheatc: cdfio.o  cdfmxlheatc.f90
	$(F90) cdfmxlheatc.f90 -o cdfmxlheatc cdfio.o $(FFLAGS)

cdfmxlheatc-full: cdfio.o  cdfmxlheatc-full.f90
	$(F90) cdfmxlheatc-full.f90 -o cdfmxlheatc-full cdfio.o $(FFLAGS)

cdfmxlsaltc: cdfio.o  cdfmxlsaltc.f90
	$(F90) cdfmxlsaltc.f90 -o cdfmxlsaltc cdfio.o $(FFLAGS)

cdfmxlhcsc: cdfio.o  eos.o cdfmxlhcsc.f90
	$(F90) cdfmxlhcsc.f90 -o cdfmxlhcsc cdfio.o eos.o $(FFLAGS)

cdficediags: cdfio.o  cdficediags.f90
	$(F90) cdficediags.f90 -o cdficediags cdfio.o $(FFLAGS)

cdfzonalmean: cdfio.o  cdfzonalmean.f90
	$(F90) cdfzonalmean.f90 -o cdfzonalmean cdfio.o $(FFLAGS) 

cdfzonalsum: cdfio.o  cdfzonalsum.f90
	$(F90) cdfzonalsum.f90 -o cdfzonalsum cdfio.o $(FFLAGS) 

cdfzonalout: cdfio.o  cdfzonalout.f90
	$(F90) cdfzonalout.f90 -o cdfzonalout cdfio.o $(FFLAGS) 

cdfhflx: cdfio.o  cdfhflx.f90
	$(F90) cdfhflx.f90 -o cdfhflx cdfio.o $(FFLAGS)

cdfwflx: cdfio.o  cdfwflx.f90
	$(F90) cdfwflx.f90 -o cdfwflx cdfio.o $(FFLAGS)

cdfbuoyflx: cdfio.o  eos.o cdfbuoyflx.f90
	$(F90) cdfbuoyflx.f90 -o cdfbuoyflx cdfio.o eos.o $(FFLAGS)

## Extracting tools, information tools
cdfprofile: cdfio.o  cdfprofile.f90
	$(F90) cdfprofile.f90  -o cdfprofile cdfio.o $(FFLAGS)

cdfwhereij: cdfio.o  cdfwhereij.f90
	$(F90) cdfwhereij.f90  -o cdfwhereij cdfio.o $(FFLAGS)

cdffindij: cdfio.o  cdffindij.f90
	$(F90) cdffindij.f90  -o cdffindij cdfio.o $(FFLAGS)

cdfweight: cdfio.o  cdfweight.f90
	$(F90) cdfweight.f90  -o cdfweight cdfio.o $(FFLAGS)

cdfcoloc: cdfio.o  cdfcoloc.f90
	$(F90) cdfcoloc.f90  -o cdfcoloc cdfio.o $(FFLAGS)

cdfcoloc2: cdfio.o  cdfcoloc2.f90
	$(F90) cdfcoloc2.f90  -o cdfcoloc2 cdfio.o $(FFLAGS)

cdfstatcoord: cdfio.o  cdfstatcoord.f90
	$(F90) cdfstatcoord.f90  -o cdfstatcoord cdfio.o $(FFLAGS)

cdfmaxmoc: cdfio.o  cdfmaxmoc.f90
	$(F90) cdfmaxmoc.f90  -o cdfmaxmoc cdfio.o $(FFLAGS)

cdfcensus: cdfio.o eos.o cdfcensus.f90
	$(F90) cdfcensus.f90  -o cdfcensus cdfio.o eos.o $(FFLAGS)

cdfzoom: cdfio.o  cdfzoom.f90
	$(F90) cdfzoom.f90  -o cdfzoom cdfio.o $(FFLAGS)

cdfmax: cdfio.o  cdfmax.f90
	$(F90) cdfmax.f90  -o cdfmax cdfio.o $(FFLAGS)

cdfmax_sp: cdfio.o  cdfmax_sp.f90
	$(F90) cdfmax_sp.f90  -o cdfmax_sp cdfio.o $(FFLAGS)

cdfprobe: cdfio.o  cdfprobe.f90
	$(F90) cdfprobe.f90  -o cdfprobe cdfio.o $(FFLAGS)

cdfclip: cdfio.o  cdfclip.f90
	$(F90) cdfclip.f90  -o cdfclip cdfio.o $(FFLAGS)

cdfsmooth: cdfio.o  cdfsmooth.f90
	$(F90) cdfsmooth.f90  -o cdfsmooth cdfio.o $(FFLAGS)

cdfpendep: cdfio.o  cdfpendep.f90
	$(F90) cdfpendep.f90  -o cdfpendep cdfio.o $(FFLAGS)

## reformating programs
cdf16bit: cdfio.o cdf16bit.f90
	$(F90) cdf16bit.f90  -o cdf16bit cdfio.o $(FFLAGS)

cdfvita: cdfio.o cdfvita.f90
	$(F90) cdfvita.f90  -o cdfvita cdfio.o $(FFLAGS)

cdftrp_bathy: cdfio.o cdftrp_bathy.f90
	$(F90) cdftrp_bathy.f90  -o cdftrp_bathy cdfio.o $(FFLAGS)

cdfconvert: cdfio.o cdfconvert.f90
	$(F90)  cdfconvert.f90  -o cdfconvert cdfio.o $(FFLAGS)

cdfflxconv: cdfio.o cdfflxconv.f90
	$(F90)   cdfflxconv.f90  -o cdfflxconv cdfio.o $(FFLAGS)

cdfsstconv: cdfio.o cdfsstconv.f90
	$(F90)   cdfsstconv.f90  -o cdfsstconv cdfio.o $(FFLAGS)

cdfstrconv: cdfio.o cdfstrconv.f90
	$(F90)   cdfstrconv.f90  -o cdfstrconv cdfio.o $(FFLAGS)

cdfbathy: cdfio.o cdfbathy.f90
	$(F90)   cdfbathy.f90  -o cdfbathy cdfio.o $(FFLAGS)

cdfcofdis: cdfio.o cdfcofdis.f90
	$(F90)    cdfcofdis.f90  -o cdfcofdis cdfio.o $(FFLAGS)

cdfcoastline: cdfio.o cdfcoastline.f90
	$(F90)    cdfcoastline.f90  -o cdfcoastline cdfio.o $(FFLAGS)

cdfvar: cdfio.o cdfvar.f90
	$(F90)   cdfvar.f90  -o cdfvar cdfio.o $(FFLAGS)

cdfcsp: cdfio.o cdfcsp.f90
	$(F90)   cdfcsp.f90  -o cdfcsp cdfio.o $(FFLAGS)

cdfpolymask: cdfio.o modpoly.o cdfpolymask.f90
	$(F90)   cdfpolymask.f90  -o cdfpolymask cdfio.o modpoly.o $(FFLAGS)


# OLD bimg/dimg stuff: use by the trpsig monitoring....
bimgmoy4: bimgmoy4.f90
	$(F90) bimgmoy4.f90  -o bimgmoy4 $(FFLAGS)

bimgcaltrans: bimgcaltrans.f90
	$(F90) bimgcaltrans.f90  -o bimgcaltrans $(FFLAGS)

coordinates2hgr: coordinates2hgr.f90
	$(F90) coordinates2hgr.f90 -o coordinates2hgr $(FFLAGS)

coordinates2zgr: coordinates2zgr.f90
	$(F90) coordinates2zgr.f90 -o coordinates2zgr $(FFLAGS)

coordinates2hgr_karine: coordinates2hgr_karine.f90
	$(F90) coordinates2hgr_karine.f90 -o coordinates2hgr_karine $(FFLAGS)

coordinates2zgr_karine: coordinates2zgr_karine.f90
	$(F90) coordinates2zgr_karine.f90 -o coordinates2zgr_karine $(FFLAGS)

## Modules

cdfio.o: cdfio.f90
	$(F90) -c  cdfio.f90 $(FFLAGS)

eos.o: eos.f90
	$(F90) -c eos.f90 $(FFLAGS)

modpoly.o: modpoly.f90
	$(F90) -c modpoly.f90 $(FFLAGS)

## Utilities
tar:
	( cd ../ ; tar cf cdftools-2.1.tar $(CDFTOOLS)/*90 $(CDFTOOLS)/Make* \
          $(CDFTOOLS)/section.dat $(CDFTOOLS)/JOBS $(CDFTOOLS)/DOC \
          $(CDFTOOLS)/macro.* )

clean:
	\rm -f *.mod *.o  *~

cleanexe: clean
	\rm -f $(EXEC)

install:
	\cp $(EXEC)  $(INSTALL)
