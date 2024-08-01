#!/bin/bash
TOP=../wt-cov-complex-solv.prmtop
INPTRAJ=../wt-cov-complex-1000ns.dcd

cpptraj <<EOF
set topname=$TOP
set inpfile=$INPTRAJ

parm \$topname
trajin \$inpfile 1 last 10
hbond HB :495,289,223,224,264 out nhb.dat avgout solute_avg.dat \
 solventacceptor :WAT@O solventdonor :WAT \
 solvout solvent_avg.dat bridgeout bridge.dat \


EOF


