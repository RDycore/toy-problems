

## Generating meshes

Run ``../share/meshes/generate_mesh_dam_break.m`` and change dx = 100, 10, 5, 2, 1.

## Run Simulations

```
mpiexec -n 1 ./ex2b -save true -output_prefix dambreak5x10 -mesh ../share/meshes/DamBreak_grid5x10.exo -initial_condition ../share/initial_conditions/DamBreak_grid5x10_wetdownstream.IC -dt 0.1 -ts_max_time 72
mpiexec -n 5 ./ex2b -save true -output_prefix dambreak50x100 -mesh ../share/meshes/DamBreak_grid50x100.exo -initial_condition ../share/initial_conditions/DamBreak_grid50x100_wetdownstream.IC -dt 0.1 -ts_max_time 72
mpiexec -n 5 ./ex2b -save true -output_prefix dambreak100x200 -mesh ../share/meshes/DamBreak_grid100x200.exo -initial_condition ../share/initial_conditions/DamBreak_grid100x200_wetdownstream.IC -dt 0.1 -ts_max_time 72
mpiexec -n 10 ./ex2b -save true -output_prefix dambreak250x500 -mesh ../share/meshes/DamBreak_grid250x500.exo -initial_condition ../share/initial_conditions/DamBreak_grid250x500_wetdownstream.IC -dt 0.1 -ts_max_time 72
mpiexec -n 40 ./ex2b -save true -output_prefix dambreak500x1000 -mesh ../share/meshes/DamBreak_grid500x1000.exo -initial_condition ../share/initial_conditions/DamBreak_grid500x1000_wetdownstream.IC -dt 0.1 -ts_max_time 72
```

## Run MMS simulations

```
mpiexec -n 1 ./ex2b_MMS -savef true -output_prefix MMS_dx1 -mesh ../share/meshes/MMS_mesh_dx1.exo -initial_condition ../share/initial_conditions/MMS_dx1.IC -mannings_n_file ../share/auxfiles/manning_n_dx1.bin -dt 0.01 -ts_max_time 5
mpiexec -n 2 ./ex2b_MMS -savef true -output_prefix MMS_dx0.5 -mesh ../share/meshes/MMS_mesh_dx0.5.exo -initial_condition ../share/initial_conditions/MMS_dx0.5.IC -mannings_n_file ../share/auxfiles/manning_n_dx0.5.bin -dt 0.01 -ts_max_time 5
mpiexec -n 2 ./ex2b_MMS -savef true -output_prefix MMS_dx0.25 -mesh ../share/meshes/MMS_mesh_dx0.25.exo -initial_condition ../share/initial_conditions/MMS_dx0.25.IC -mannings_n_file ../share/auxfiles/manning_n_dx0.25.bin -dt 0.01 -ts_max_time 5
mpiexec -n 2 ./ex2b_MMS -savef true -output_prefix MMS_dx0.1 -mesh ../share/meshes/MMS_mesh_dx0.1.exo -initial_condition ../share/initial_conditions/MMS_dx0.1.IC -mannings_n_file ../share/auxfiles/manning_n_dx0.1.bin -dt 0.01 -ts_max_time 5
mpiexec -n 4 ./ex2b_MMS -savef true -output_prefix MMS_dx0.05 -mesh ../share/meshes/MMS_mesh_dx0.05.exo -initial_condition ../share/initial_conditions/MMS_dx0.05.IC -mannings_n_file ../share/auxfiles/manning_n_dx0.05.bin -dt 0.01 -ts_max_time 5
```
