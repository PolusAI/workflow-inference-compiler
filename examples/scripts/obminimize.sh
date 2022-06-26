# -cg is the default minimization algorithm, but it falsely claims to reach
# convergence on the second timestep. (same with -newton)
# Use steepest descent because it's the only one that works!
obminimize -sd -o mol2 $1 > temp.mol2 # -xu is ignored here
obabel temp.mol2 -o mol2 -xu