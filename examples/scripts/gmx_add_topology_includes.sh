#!/bin/bash -e
# NOTE: These include statements MUST be in this order. In particular, we
# cannot simply use cat $1 because we need to interleave the amber99sb-ildn.ff
# includes around ligand_GMX.itp, and because we need to remove the [ defaults ]
# directive (which is included in forcefield.itp)
echo '#include "amber99sb-ildn.ff/forcefield.itp"'

echo '#include "ligand_GMX.itp"'

echo '#include "amber99sb-ildn.ff/spce.itp"'
echo '#include "amber99sb-ildn.ff/ions.itp"'

echo "[ system ]"
echo " ligand "
echo ""
echo "[ molecules ]"
echo "; Compound        nmols"
echo " ligand           1     "

#cat $1
