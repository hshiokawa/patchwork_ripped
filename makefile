
CC=mpicc -cc=icc

FLAGS = -lhdf5_hl -lhdf5 
I_LOCS = -I/usr/local/encap/hdf5-1.8.4/include -c
L_LOCS = -L/usr/local/encap/hdf5-1.8.4/lib

CC_COMPILE  = $(CC) $(FLAGS) $(I_LOCS) 
CC_LOAD     = $(CC) $(FLAGS) $(L_LOCS)



#############################################################

.c.o:
	$(CC_COMPILE) $*.c

EXE  = harm3d_patch

all: $(EXE)



#############################################################

SRCS = \
patchwork_gl.c patchwork_lg.c patchwork_setup.c patchwork_move.c \
misc.c coord.c harm_mpi.c init.c decomp.c interp.c dump_ascii.c \
dump_hdf.c diag.c restart.c metric.c fixup.c bounds.c step_ch.c \
recon.c phys.c utoprim_2d_fast.c u2p_util.c utoprim_1d.c \
transform.c lu.c conn_func.c dxpdxp_calc.c utoprim_1d_ee.c \
utoprim_1d_ee2.c timestep.c


BASE = $(basename $(SRCS) )
OBJS = $(addsuffix .o, $(BASE) )

INCS = \
decs.h patchwork.h defs.h patchwork_defs.h harm_mpi.h metric.h \
recon.h u2p_util.h u2p_defs.h



#############################################################

defs.h: decs.h 
	egrep "^extern|^#if|^#elif|^#else|^#endif" decs.h | sed s/^extern//  > defs.h

patchwork_defs.h: patchwork.h 
	egrep "^extern|^#if|^#elif|^#else|^#endif" patchwork.h | sed s/^extern//  > patchwork_defs.h

$(OBJS): $(INCS) makefile

main.o: $(INCS) 



#############################################################

$(EXE): main.o $(OBJS)
	$(CC_LOAD) $(OBJS) main.o -o $(EXE)



#############################################################

clean:
	/bin/rm -f *.o
	/bin/rm -f defs.h 
	/bin/rm -f patchwork_defs.h 
	/bin/rm -f $(EXE)
