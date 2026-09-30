K_SRC=TDMA_solver
HOST_SRC=adi_app.cpp xcl2.cpp golden_fun.cpp
CONFIG=TDMA.ini
APP_ARGS=-nx=32 -ny=32 -nz=32 -iters=1


include ../../utility/Makefile.mk

# The kernel includes "pre_proc.cpp" from this directory: the sw_emu link recompiles the
# kernel from a copy under build/, so the directory must be an explicit include path.
KCC_INCDIR += -I.
