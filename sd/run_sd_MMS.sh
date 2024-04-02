source /qfs/people/feng779/RDycore/petsc-rdycore-intel.sh
suffixes=("1" "0.5" "0.25" "0.1" "0.05")
for suffix in "${suffixes[@]}"
do
  mpiexec -n 2 \
  ./ex2d_MMS \
  -savef true \
  -output_prefix MMS_dx${suffix} \
  -mesh ../share/meshes/ex2b_MMS_mesh_dx${suffix}.exo \
  -sed 1 \
  -dt 0.01 \
  -ts_max_time 5.0
done

