#This has been tested in Bash, though I'm rather confident that it will go on zsh as well.

#xtb things
export XTBHOME=$BTROOT/xtb-6.6.0
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=2,1 #I'm assuming 2 cores. Feel free to change this
export MKL_NUM_THREADS=2
ulimit -s unlimited
source $XTBHOME/share/xtb/config_env.bash

#Not strictly needed, but nice to have
export PATH=$BTROOT:$PATH

#xdrlib
export LD_LIBRARY_PATH=$BTROOT/xdrlib/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$BTROOT/xdrlib/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$BTROOT/xdrlib/include:$CPLUS_INCLUDE_PATH
