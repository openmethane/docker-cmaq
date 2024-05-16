#   Makefile generated using program [bldmake]
#
#   Generation date [06/30/2022  18:02:08]
#   Configuration file [cfg.bwd_CH4only.bld]
#   Using GIT repository
#

MODEL = ADJOINT_BWD

FC = mpif90
CC = mpicc

ROOT = $(PWD)/..
IOAPI_ROOT = $(PWD)/../../ioapi-3.1
REPOROOT = $(PWD)/../modified_CCTM
CMAQ_ROOT = "/opt/cmaq/CMAQ-5.0.2"

f_FLAGS    = -extend_source 132 -qopt-report0 -nodefines -warn nousage -stand none -traceback -CB -c -g -O2 -module ../BLD_bwd_CH4only/MOD_DIR -I.
f_FLAGS    = -ffixed-line-length-none -fd-lines-as-comments -cpp -O2 -J../BLD_fwd_CH4only/MOD_DIR -I. -I$(CMAQ_ROOT)/lib/$(ARCH)/gcc/se_snl -I$(CMAQ_ROOT)/lib/$(ARCH)/gcc/pario
F_FLAGS    = $(f_FLAGS)
f90_FLAGS  = 
F90_FLAGS  = $(f90_FLAGS)
C_FLAGS    =

LINKER     = $(FC)
LINK_FLAGS = -qopenmp


CPP = $(FC)
CPP_FLAGS = \
  -Dparallel \
  -Dcluster \
  -DINTERPB=PINTERPB \
  -DM3ERR=PM3ERR \
  -DM3EXIT=PM3EXIT \
  -DM3WARN=PM3WARN \
  -DSHUT3=PSHUT3 \
  -DWRITE3=PWRITE3 \
  -DSUBST_MODULES=SE_MODULES \
  -DSUBST_BARRIER=SE_BARRIER \
  -DSUBST_GLOBAL_MAX=SE_GLOBAL_MAX \
  -DSUBST_GLOBAL_MIN=SE_GLOBAL_MIN \
  -DSUBST_GLOBAL_MIN_DATA=SE_GLOBAL_MIN_DATA \
  -DSUBST_GLOBAL_TO_LOCAL_COORD=SE_GLOBAL_TO_LOCAL_COORD \
  -DSUBST_GLOBAL_SUM=SE_GLOBAL_SUM \
  -DSUBST_GLOBAL_LOGICAL=SE_GLOBAL_LOGICAL \
  -DSUBST_LOOP_INDEX=SE_LOOP_INDEX \
  -DSUBST_SUBGRID_INDEX=SE_SUBGRID_INDEX \
  -DSUBST_HI_LO_BND_PE=SE_HI_LO_BND_PE \
  -DSUBST_SUM_CHK=SE_SUM_CHK \
  -DSUBST_INIT_ARRAY=SE_INIT_ARRAY \
  -DSUBST_COMM=SE_COMM \
  -DSUBST_MY_REGION=SE_MY_REGION \
  -DSUBST_SLICE=SE_SLICE \
  -DSUBST_GATHER=SE_GATHER \
  -DSUBST_DATA_COPY=SE_DATA_COPY \
  -DSUBST_IN_SYN=SE_IN_SYN

SE_SNL =  -L$(CMAQ_ROOT)/lib/$(ARCH)/gcc/se_snl -lse_snl
PARIO =  -L$(CMAQ_ROOT)/lib/$(ARCH)/gcc/pario -lpario
MPI = -L/usr/lib/aarch64-linux-gnu/openmpi/lib -lmpi
MPI_MPIFH =  -L/usr/lib/aarch64-linux-gnu/openmpi/lib -lmpi_mpifh
IOAPI =  -L$(IOAPI_ROOT)/$BIN -lioapi
NETCDF =  -lnetcdf
NETCDFF =  -lnetcdff

LIBRARIES = $(SE_SNL) $(PARIO) $(MPI) $(MPI_MPIFH) $(IOAPI) $(NETCDF) $(NETCDFF)

BASE_INC = $(ROOT)/ICL/fixed
MECH_INC = $(ROOT)/ICL/mech/cb05cl_ae5_aq_CH4only
PA_INC = $(ROOT)/ICL/procan/pa_noop

INCLUDES =  \
  -DSUBST_PE_COMM=\"$(BASE_INC)/par/PE_COMM.EXT\" \
  -DSUBST_CONST=\"$(BASE_INC)/const/CONST.EXT\" \
  -DSUBST_FILES_ID=\"$(BASE_INC)/filenames/FILES_CTM.EXT\" \
  -DSUBST_EMPR_VD=\"$(BASE_INC)/emctrl/EMISPRM.vdif.EXT\" \
  -DSUBST_EMPR_CH=\"$(BASE_INC)/emctrl/EMISPRM.chem.EXT\" \
  -DSUBST_IOPARMS=\"$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT\" \
  -DSUBST_IOFDESC=\"$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT\" \
  -DSUBST_IODECL=\"$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT\" \
  -DSUBST_RXCMMN=\"$(MECH_INC)/RXCM.EXT\" \
  -DSUBST_RXDATA=\"$(MECH_INC)/RXDT.EXT\" \
  -DSUBST_GC_SPC=\"$(MECH_INC)/GC_SPC.EXT\" \
  -DSUBST_GC_EMIS=\"$(MECH_INC)/GC_EMIS.EXT\" \
  -DSUBST_GC_ICBC=\"$(MECH_INC)/GC_ICBC.EXT\" \
  -DSUBST_GC_DIFF=\"$(MECH_INC)/GC_DIFF.EXT\" \
  -DSUBST_GC_DDEP=\"$(MECH_INC)/GC_DDEP.EXT\" \
  -DSUBST_GC_DEPV=\"$(MECH_INC)/GC_DEPV.EXT\" \
  -DSUBST_GC_ADV=\"$(MECH_INC)/GC_ADV.EXT\" \
  -DSUBST_GC_CONC=\"$(MECH_INC)/GC_CONC.EXT\" \
  -DSUBST_GC_G2AE=\"$(MECH_INC)/GC_G2AE.EXT\" \
  -DSUBST_GC_G2AQ=\"$(MECH_INC)/GC_G2AQ.EXT\" \
  -DSUBST_GC_SCAV=\"$(MECH_INC)/GC_SCAV.EXT\" \
  -DSUBST_GC_WDEP=\"$(MECH_INC)/GC_WDEP.EXT\" \
  -DSUBST_AE_SPC=\"$(MECH_INC)/AE_SPC.EXT\" \
  -DSUBST_AE_EMIS=\"$(MECH_INC)/AE_EMIS.EXT\" \
  -DSUBST_AE_ICBC=\"$(MECH_INC)/AE_ICBC.EXT\" \
  -DSUBST_AE_DIFF=\"$(MECH_INC)/AE_DIFF.EXT\" \
  -DSUBST_AE_DDEP=\"$(MECH_INC)/AE_DDEP.EXT\" \
  -DSUBST_AE_DEPV=\"$(MECH_INC)/AE_DEPV.EXT\" \
  -DSUBST_AE_ADV=\"$(MECH_INC)/AE_ADV.EXT\" \
  -DSUBST_AE_CONC=\"$(MECH_INC)/AE_CONC.EXT\" \
  -DSUBST_AE_A2AQ=\"$(MECH_INC)/AE_A2AQ.EXT\" \
  -DSUBST_AE_SCAV=\"$(MECH_INC)/AE_SCAV.EXT\" \
  -DSUBST_AE_WDEP=\"$(MECH_INC)/AE_WDEP.EXT\" \
  -DSUBST_NR_SPC=\"$(MECH_INC)/NR_SPC.EXT\" \
  -DSUBST_NR_EMIS=\"$(MECH_INC)/NR_EMIS.EXT\" \
  -DSUBST_NR_ICBC=\"$(MECH_INC)/NR_ICBC.EXT\" \
  -DSUBST_NR_DIFF=\"$(MECH_INC)/NR_DIFF.EXT\" \
  -DSUBST_NR_DDEP=\"$(MECH_INC)/NR_DDEP.EXT\" \
  -DSUBST_NR_DEPV=\"$(MECH_INC)/NR_DEPV.EXT\" \
  -DSUBST_NR_ADV=\"$(MECH_INC)/NR_ADV.EXT\" \
  -DSUBST_NR_N2AE=\"$(MECH_INC)/NR_N2AE.EXT\" \
  -DSUBST_NR_CONC=\"$(MECH_INC)/NR_CONC.EXT\" \
  -DSUBST_NR_N2AQ=\"$(MECH_INC)/NR_N2AQ.EXT\" \
  -DSUBST_NR_SCAV=\"$(MECH_INC)/NR_SCAV.EXT\" \
  -DSUBST_NR_WDEP=\"$(MECH_INC)/NR_WDEP.EXT\" \
  -DSUBST_TR_SPC=\"$(ROOT)/ICL/mech/trac0/TR_SPC.EXT\" \
  -DSUBST_TR_EMIS=\"$(ROOT)/ICL/mech/trac0/TR_EMIS.EXT\" \
  -DSUBST_TR_ICBC=\"$(ROOT)/ICL/mech/trac0/TR_ICBC.EXT\" \
  -DSUBST_TR_DIFF=\"$(ROOT)/ICL/mech/trac0/TR_DIFF.EXT\" \
  -DSUBST_TR_DDEP=\"$(ROOT)/ICL/mech/trac0/TR_DDEP.EXT\" \
  -DSUBST_TR_DEPV=\"$(ROOT)/ICL/mech/trac0/TR_DEPV.EXT\" \
  -DSUBST_TR_ADV=\"$(ROOT)/ICL/mech/trac0/TR_ADV.EXT\" \
  -DSUBST_TR_T2AQ=\"$(ROOT)/ICL/mech/trac0/TR_T2AQ.EXT\" \
  -DSUBST_TR_SCAV=\"$(ROOT)/ICL/mech/trac0/TR_SCAV.EXT\" \
  -DSUBST_TR_WDEP=\"$(ROOT)/ICL/mech/trac0/TR_WDEP.EXT\" \
  -DSUBST_PACTL_ID=\"$(PA_INC)/PA_CTL.EXT\" \
  -DSUBST_PACMN_ID=\"$(PA_INC)/PA_CMN.EXT\" \
  -DSUBST_PADAT_ID=\"$(PA_INC)/PA_DAT.EXT\" \
  -DSUBST_MPICH=\"$(OPENMPI)\"

GLOBAL_MODULES = \
  CGRID_SPCS.o \
  HGRD_DEFN.o \
  VGRD_DEFN.o \
  GRID_CONF.o \
  AERO_EMIS.o \
  BIOG_EMIS.o \
  adj_array_mod.o \
  BEIS_DEFN.o \
  PMEM_DEFN.o \
  UDTYPES.o \
  STK_EMIS.o \
  STK_PRMS.o \
  PTMAP.o \
  PTBILIN.o \
  PTMET.o \
  PT3D_DEFN.o \
  EMIS_DEFN.o \
  VDIFF_MAP.o \
  PCGRID_DEFN.o \
  CHK_FILES_DEFS.o \
  PAGRD_DEFN.o \
  CHK_FILES_IMPL.o \
  ADJOINT_FILES.o \
  CHK_FILES.o

DRIVER = \
  advstep.o \
  driver_bwd.o \
  hcontvel.o \
  hveloc.o \
  readconf.o \
  sciproc_bwd.o \
  syncstep.o

PAR = \
  distr_env.o \
  par_init.o \
  par_term.o

INIT = \
  define_receptor.o \
  diffmsg.o \
  flcheck.o \
  grdcheck.o \
  init_lgrid.o \
  setup_logdev.o

COUPLE = \
  couple.o \
  ppmv_msmr.o

HADV = \
  hadv_bwd.o \
  hppm.o \
  rdbcon_rj.o \
  xadv_bwd.o \
  xadv_rj_fwd.o \
  yadv_bwd.o \
  yadv_rj_fwd.o \
  zfdbc.o

VADV = \
  vppm_bwd.o \
  zadvyppm_bwd.o

HDIFF = \
  deform.o \
  hcdiff3d.o \
  hdiff_bwd.o \
  rho_j.o

VDIFF = \
  DEPVVARS.o \
  DEPV_DEFN.o \
  beis3.o \
  cgrid_depv.o \
  checkmem.o \
  chkgrid.o \
  conv_cgrid_bwd.o \
  conv_em_bwd.o \
  czangle.o \
  delta_zs.o \
  eddyx.o \
  fire_plmris.o \
  gas_depv_map.o \
  getfline.o \
  getparb.o \
  hrno.o \
  m3dry.o \
  matrix_bwd.o \
  opddep.o \
  opdepv_diag.o \
  opemis.o \
  openlayout.o \
  opmet.o \
  oppt3d_diag.o \
  parsline.o \
  plmris.o \
  plsprd.o \
  preplm.o \
  rddepv.o \
  rdmet.o \
  tmpbeis.o \
  tri_bwd.o \
  ungridb2.o \
  vdiffacm2_bwd.o \
  wrdaymsg.o \
  write3_distr.o

PHOT = \
  phot_noop.o

CHEM = \
  chem_bwd_noop.o

AERO = \
  aero_bwd_noop.o

AERO_DEPV = \
  aero_depv_noop.o

CLOUD = \
  cloud_noop.o \
  hlconst.o

PROCAN = \
  pa_init.o \
  pa_irr.o \
  pa_irr_ctl.o \
  pa_mkhdr.o \
  pa_output.o \
  pa_update.o

UTIL = \
  adBuffer.o \
  adStack.o \
  cksummer.o \
  findex.o \
  get_envlist.o \
  lstepf.o \
  subhdomain.o \
  subhfile.o

OBJS = \
  $(GLOBAL_MODULES) \
  $(DRIVER) \
  $(PAR) \
  $(INIT) \
  $(COUPLE) \
  $(HADV) \
  $(VADV) \
  $(HDIFF) \
  $(VDIFF) \
  $(PHOT) \
  $(CHEM) \
  $(AERO) \
  $(AERO_DEPV) \
  $(CLOUD) \
  $(PROCAN) \
  $(UTIL)

.SUFFIXES: .F .f .c .F90 .f90

$(MODEL): $(OBJS)
	$(LINKER) $(LINK_FLAGS) $(OBJS) $(LIBRARIES) -o $@

.F.o:
	$(FC) -c $(F_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f.o:
	$(FC) -c $(f_FLAGS) $<

.F90.o:
	$(FC) -c $(F90_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f90.o:
	$(FC) -c $(f90_FLAGS) $<

.c.o:
	$(CC) -c $(C_FLAGS) $<

clean:
	rm -f $(OBJS) **/*.mod


# dependencies

ADJOINT_FILES.o:	GRID_CONF.F CHK_FILES_DEFS.F EMIS_DEFN.F VDIFF_MAP.F \
		CGRID_SPCS.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_DIFF.EXT $(MECH_INC)/AE_SPC.EXT \
		$(MECH_INC)/AE_DIFF.EXT $(MECH_INC)/NR_SPC.EXT $(MECH_INC)/NR_DIFF.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DIFF.EXT
CGRID_SPCS.o:	$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/NR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT
GRID_CONF.o:	HGRD_DEFN.F VGRD_DEFN.F
HGRD_DEFN.o:	$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
PCGRID_DEFN.o:	GRID_CONF.F CGRID_SPCS.F EMIS_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_DIFF.EXT $(MECH_INC)/AE_DIFF.EXT $(MECH_INC)/NR_DIFF.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DIFF.EXT
VGRD_DEFN.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
advstep.o:	GRID_CONF.F $(BASE_INC)/const/CONST.EXT $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
driver_bwd.o:	PCGRID_DEFN.F VGRD_DEFN.F ADJOINT_FILES.F adj_array_mod.f \
		CHK_FILES_DEFS.F CHK_FILES.F EMIS_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/RXCM.EXT $(MECH_INC)/RXDT.EXT $(MECH_INC)/GC_SPC.EXT \
		$(MECH_INC)/AE_SPC.EXT $(MECH_INC)/NR_SPC.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT
hcontvel.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
hveloc.o:	GRID_CONF.F $(BASE_INC)/par/PE_COMM.EXT $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
sciproc_bwd.o:	adj_array_mod.f $(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
syncstep.o:	GRID_CONF.F $(BASE_INC)/const/CONST.EXT $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
par_init.o:	GRID_CONF.F PAGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(PA_INC)/PA_CTL.EXT $(OPENMPI)
par_term.o:	$(OPENMPI)
define_receptor.o:	GRID_CONF.F CGRID_SPCS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/NR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT
diffmsg.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
flcheck.o:	$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
grdcheck.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
init_lgrid.o:	ADJOINT_FILES.F adj_array_mod.f $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
setup_logdev.o:	$(OPENMPI)
couple.o:	GRID_CONF.F CGRID_SPCS.F CHK_FILES_DEFS.F CHK_FILES.F \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/NR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT
ppmv_msmr.o:	GRID_CONF.F $(BASE_INC)/const/CONST.EXT $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
hadv_bwd.o:	GRID_CONF.F CHK_FILES_DEFS.F CHK_FILES.F $(BASE_INC)/par/PE_COMM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(MECH_INC)/GC_ADV.EXT $(MECH_INC)/AE_ADV.EXT $(MECH_INC)/NR_ADV.EXT \
		$(ROOT)/ICL/mech/trac0/TR_ADV.EXT
hppm.o:	HGRD_DEFN.F
rdbcon_rj.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
xadv_bwd.o:	GRID_CONF.F CGRID_SPCS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_ADV.EXT $(MECH_INC)/AE_ADV.EXT \
		$(MECH_INC)/NR_ADV.EXT $(ROOT)/ICL/mech/trac0/TR_ADV.EXT
xadv_rj_fwd.o:	HGRD_DEFN.F
yadv_bwd.o:	GRID_CONF.F CGRID_SPCS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_ADV.EXT $(MECH_INC)/AE_ADV.EXT \
		$(MECH_INC)/NR_ADV.EXT $(ROOT)/ICL/mech/trac0/TR_ADV.EXT
yadv_rj_fwd.o:	HGRD_DEFN.F
vppm_bwd.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT
zadvyppm_bwd.o:	GRID_CONF.F CGRID_SPCS.F CHK_FILES_DEFS.F CHK_FILES.F \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_ADV.EXT $(MECH_INC)/AE_SPC.EXT \
		$(MECH_INC)/AE_ADV.EXT $(MECH_INC)/NR_SPC.EXT $(MECH_INC)/NR_ADV.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_ADV.EXT
deform.o:	GRID_CONF.F $(BASE_INC)/par/PE_COMM.EXT $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
hcdiff3d.o:	GRID_CONF.F $(BASE_INC)/par/PE_COMM.EXT $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
hdiff_bwd.o:	GRID_CONF.F CGRID_SPCS.F $(BASE_INC)/const/CONST.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_DIFF.EXT $(MECH_INC)/AE_DIFF.EXT $(MECH_INC)/NR_DIFF.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DIFF.EXT
rho_j.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
AERO_EMIS.o:	GRID_CONF.F $(MECH_INC)/AE_EMIS.EXT
BEIS_DEFN.o:	HGRD_DEFN.F BIOG_EMIS.F $(BASE_INC)/const/CONST.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/NR_EMIS.EXT $(ROOT)/ICL/mech/trac0/TR_EMIS.EXT
BIOG_EMIS.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
DEPVVARS.o:	HGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT
DEPV_DEFN.o:	HGRD_DEFN.F DEPVVARS.F GRID_CONF.F CGRID_SPCS.F $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(MECH_INC)/RXCM.EXT $(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_DEPV.EXT \
		$(MECH_INC)/AE_SPC.EXT $(MECH_INC)/AE_DEPV.EXT $(MECH_INC)/NR_SPC.EXT \
		$(MECH_INC)/NR_DEPV.EXT $(MECH_INC)/NR_N2AE.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DEPV.EXT
EMIS_DEFN.o:	GRID_CONF.F BEIS_DEFN.F PT3D_DEFN.F PMEM_DEFN.F AERO_EMIS.F \
		BIOG_EMIS.F PTMAP.F $(BASE_INC)/const/CONST.EXT $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/AE_EMIS.EXT \
		$(MECH_INC)/NR_EMIS.EXT $(ROOT)/ICL/mech/trac0/TR_EMIS.EXT
PMEM_DEFN.o:	HGRD_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT
PT3D_DEFN.o:	UDTYPES.F GRID_CONF.F STK_PRMS.F STK_EMIS.F PMEM_DEFN.F \
		PTMET.F PTBILIN.F PTMAP.F VGRD_DEFN.F $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/AE_EMIS.EXT $(MECH_INC)/NR_EMIS.EXT \
		$(ROOT)/ICL/mech/trac0/TR_EMIS.EXT
PTBILIN.o:	UDTYPES.F HGRD_DEFN.F VGRD_DEFN.F STK_PRMS.F
PTMAP.o:	UDTYPES.F STK_EMIS.F PMEM_DEFN.F $(MECH_INC)/GC_EMIS.EXT \
		$(MECH_INC)/AE_EMIS.EXT $(MECH_INC)/NR_EMIS.EXT $(ROOT)/ICL/mech/trac0/TR_EMIS.EXT
PTMET.o:	UDTYPES.F STK_PRMS.F PTBILIN.F VGRD_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
STK_EMIS.o:	UDTYPES.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
STK_PRMS.o:	UDTYPES.F HGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
VDIFF_MAP.o:	CGRID_SPCS.F AERO_EMIS.F $(BASE_INC)/emctrl/EMISPRM.vdif.EXT \
		$(BASE_INC)/emctrl/EMISPRM.chem.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/GC_DIFF.EXT \
		$(MECH_INC)/GC_DDEP.EXT $(MECH_INC)/GC_DEPV.EXT $(MECH_INC)/AE_SPC.EXT \
		$(MECH_INC)/AE_DIFF.EXT $(MECH_INC)/AE_DDEP.EXT $(MECH_INC)/AE_DEPV.EXT \
		$(MECH_INC)/NR_SPC.EXT $(MECH_INC)/NR_EMIS.EXT $(MECH_INC)/NR_DIFF.EXT \
		$(MECH_INC)/NR_DDEP.EXT $(MECH_INC)/NR_DEPV.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_EMIS.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DIFF.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DDEP.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DEPV.EXT
beis3.o:	BIOG_EMIS.F
cgrid_depv.o:	HGRD_DEFN.F $(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_DEPV.EXT \
		$(MECH_INC)/AE_SPC.EXT $(MECH_INC)/AE_DEPV.EXT $(MECH_INC)/NR_SPC.EXT \
		$(MECH_INC)/NR_DEPV.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DEPV.EXT
chkgrid.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT
conv_cgrid_bwd.o:	GRID_CONF.F CGRID_SPCS.F $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/NR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT
conv_em_bwd.o:	GRID_CONF.F CGRID_SPCS.F EMIS_DEFN.F $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_DIFF.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/AE_DIFF.EXT \
		$(MECH_INC)/NR_DIFF.EXT $(ROOT)/ICL/mech/trac0/TR_DIFF.EXT
czangle.o:	$(BASE_INC)/const/CONST.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
eddyx.o:	GRID_CONF.F $(BASE_INC)/par/PE_COMM.EXT $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
fire_plmris.o:	$(BASE_INC)/const/CONST.EXT
gas_depv_map.o:	DEPVVARS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_DEPV.EXT $(MECH_INC)/NR_SPC.EXT \
		$(MECH_INC)/NR_DEPV.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DEPV.EXT
hrno.o:	HGRD_DEFN.F BIOG_EMIS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT
m3dry.o:	HGRD_DEFN.F DEPVVARS.F $(BASE_INC)/const/CONST.EXT $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
matrix_bwd.o:	VGRD_DEFN.F $(MECH_INC)/GC_DIFF.EXT $(MECH_INC)/AE_DIFF.EXT \
		$(MECH_INC)/NR_DIFF.EXT $(ROOT)/ICL/mech/trac0/TR_DIFF.EXT
opddep.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_DDEP.EXT $(MECH_INC)/AE_DDEP.EXT $(MECH_INC)/NR_DDEP.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DDEP.EXT
opdepv_diag.o:	HGRD_DEFN.F DEPVVARS.F $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
opemis.o:	VGRD_DEFN.F EMIS_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/NR_EMIS.EXT $(ROOT)/ICL/mech/trac0/TR_EMIS.EXT
openlayout.o:	HGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
opmet.o:	$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
oppt3d_diag.o:	GRID_CONF.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
plmris.o:	$(BASE_INC)/const/CONST.EXT
rddepv.o:	HGRD_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/RXCM.EXT $(MECH_INC)/GC_DEPV.EXT $(MECH_INC)/AE_DEPV.EXT \
		$(MECH_INC)/NR_DEPV.EXT $(ROOT)/ICL/mech/trac0/TR_DEPV.EXT
rdmet.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
tmpbeis.o:	HGRD_DEFN.F BIOG_EMIS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT
tri_bwd.o:	VGRD_DEFN.F $(MECH_INC)/GC_DIFF.EXT $(MECH_INC)/AE_DIFF.EXT \
		$(MECH_INC)/NR_DIFF.EXT $(ROOT)/ICL/mech/trac0/TR_DIFF.EXT
ungridb2.o:	HGRD_DEFN.F PTBILIN.F
vdiffacm2_bwd.o:	PCGRID_DEFN.F EMIS_DEFN.F DEPV_DEFN.F VGRD_DEFN.F \
		VDIFF_MAP.F CHK_FILES_DEFS.F CHK_FILES.F adj_array_mod.f ADJOINT_FILES.F \
		$(BASE_INC)/const/CONST.EXT $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT
phot_noop.o:	HGRD_DEFN.F VGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/RXCM.EXT
chem_bwd_noop.o:	$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
aero_bwd_noop.o:	$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
aero_depv_noop.o:	HGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/AE_DEPV.EXT
cloud_noop.o:	$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
hlconst.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
PAGRD_DEFN.o:	GRID_CONF.F $(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(PA_INC)/PA_CTL.EXT $(PA_INC)/PA_CMN.EXT
pa_init.o:	GRID_CONF.F PAGRD_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(PA_INC)/PA_CTL.EXT $(PA_INC)/PA_CMN.EXT $(PA_INC)/PA_DAT.EXT
pa_irr.o:	GRID_CONF.F $(MECH_INC)/RXCM.EXT $(MECH_INC)/GC_SPC.EXT \
		$(PA_INC)/PA_CMN.EXT
pa_irr_ctl.o:	GRID_CONF.F PAGRD_DEFN.F $(MECH_INC)/GC_SPC.EXT \
		$(PA_INC)/PA_CMN.EXT
pa_mkhdr.o:	GRID_CONF.F PAGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/AE_SPC.EXT $(PA_INC)/PA_CMN.EXT
pa_output.o:	GRID_CONF.F PAGRD_DEFN.F $(BASE_INC)/filenames/FILES_CTM.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(PA_INC)/PA_CTL.EXT $(PA_INC)/PA_CMN.EXT
pa_update.o:	PAGRD_DEFN.F AERO_EMIS.F $(BASE_INC)/const/CONST.EXT \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(BASE_INC)/emctrl/EMISPRM.vdif.EXT \
		$(BASE_INC)/emctrl/EMISPRM.chem.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/GC_DDEP.EXT \
		$(MECH_INC)/GC_DEPV.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/AE_DDEP.EXT \
		$(MECH_INC)/AE_DEPV.EXT $(MECH_INC)/NR_SPC.EXT $(MECH_INC)/NR_EMIS.EXT \
		$(MECH_INC)/NR_DDEP.EXT $(MECH_INC)/NR_DEPV.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_EMIS.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DDEP.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DEPV.EXT \
		$(PA_INC)/PA_CTL.EXT $(PA_INC)/PA_CMN.EXT
cksummer.o:	GRID_CONF.F CGRID_SPCS.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/AE_SPC.EXT $(MECH_INC)/NR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_SPC.EXT
lstepf.o:	$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
subhdomain.o:	$(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
subhfile.o:	GRID_CONF.F $(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
CHK_FILES.o:	GRID_CONF.F CHK_FILES_DEFS.F CHK_FILES_IMPL.F adj_array_mod.f \
		$(BASE_INC)/filenames/FILES_CTM.EXT $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT \
		$(MECH_INC)/GC_SPC.EXT $(MECH_INC)/GC_EMIS.EXT $(MECH_INC)/GC_DIFF.EXT \
		$(MECH_INC)/GC_ADV.EXT $(MECH_INC)/GC_CONC.EXT $(MECH_INC)/AE_SPC.EXT \
		$(MECH_INC)/AE_EMIS.EXT $(MECH_INC)/AE_DIFF.EXT $(MECH_INC)/AE_ADV.EXT \
		$(MECH_INC)/AE_CONC.EXT $(MECH_INC)/NR_SPC.EXT $(MECH_INC)/NR_EMIS.EXT \
		$(MECH_INC)/NR_DIFF.EXT $(MECH_INC)/NR_ADV.EXT $(ROOT)/ICL/mech/trac0/TR_SPC.EXT \
		$(ROOT)/ICL/mech/trac0/TR_EMIS.EXT \
		$(ROOT)/ICL/mech/trac0/TR_DIFF.EXT \
		$(ROOT)/ICL/mech/trac0/TR_ADV.EXT
CHK_FILES_IMPL.o:	CHK_FILES_DEFS.F HGRD_DEFN.F VGRD_DEFN.F $(IOAPI_ROOT)/ioapi/fixed_src/PARMS3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/FDESC3.EXT \
		$(IOAPI_ROOT)/ioapi/fixed_src/IODECL3.EXT
