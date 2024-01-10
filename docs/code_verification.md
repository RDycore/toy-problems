

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
suffixes=("1" "0.5" "0.25" "0.1" "0.05")
for suffix in $suffixes
do
  ./ex2b_MMS \
  -savef true \
  -output_prefix MMS_dx${suffix} \
  -mesh ../share/meshes/ex2b_MMS_mesh_dx${suffix}.exo \
  -dt 0.01 \
  -ts_max_time 5.0 \
  -ts_exact_final_time matchstep
done
```

## Run validation cases

C property (immersed bump)
```
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_c_property_immersed -mesh ../share/meshes/ex2b_c_property_immersed.exo -initial_condition ../share/initial_conditions/ex2b_c_property_immersed.IC -mannings_n 0 -dt 0.0001 -ts_max_time 20 -caseid 1.1
```
C property (emerged bump)
```
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_c_property_emerged -mesh ../share/meshes/ex2b_c_property_emerged.exo -initial_condition ../share/initial_conditions/ex2b_c_property_emerged.IC -mannings_n 0 -dt 0.0001 -ts_max_time 20 -caseid 1.2
```

Well balancing (subcritical to supercritical)
```
mpiexec -n 20 ./ex2b_validation -savef true -output_prefix ex2b_well_balancing_sub_to_super -mesh ../share/meshes/ex2b_well_balancing_sub_to_super.exo -initial_condition ../share/initial_conditions/ex2b_well_balancing_sub_to_super.IC -mannings_n 0.0218 -dt 0.01 -ts_max_time 3000 -caseid 2.1
```
Well balancing (supercritical to subcritical)
```
mpiexec -n 20 ./ex2b_validation -savef true -output_prefix ex2b_well_balancing_super_to_sub -mesh ../share/meshes/ex2b_well_balancing_super_to_sub.exo -initial_condition ../share/initial_conditions/ex2b_well_balancing_super_to_sub.IC -mannings_n 0.0218 -dt 0.01 -ts_max_time 3000 -caseid 2.2
```

Dam break (wet domain)
```
mpiexec -n 2 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_wet_domain_dx0.1 -mesh ../share/meshes/ex2b_dam_break_wet_domain_dx0.1.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_wet_domain_dx0.1.IC -mannings_n 0 -dt 0.001 -ts_max_time 6 -caseid 3.1
mpiexec -n 4 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_wet_domain_dx0.01 -mesh ../share/meshes/ex2b_dam_break_wet_domain_dx0.01.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_wet_domain_dx0.01.IC -mannings_n 0 -dt 0.001 -ts_max_time 6 -caseid 3.1
mpiexec -n 8 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_wet_domain_dx0.001 -mesh ../share/meshes/ex2b_dam_break_wet_domain_dx0.001.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_wet_domain_dx0.001.IC -mannings_n 0 -dt 0.001 -ts_max_time 6 -caseid 3.1
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_wet_domain_dx0.0001 -mesh ../share/meshes/ex2b_dam_break_wet_domain_dx0.0001.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_wet_domain_dx0.0001.IC -mannings_n 0 -dt 0.0001 -ts_max_time 6 -caseid 3.1
```
Dam break (dry domain)
```
mpiexec -n 2 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_dry_domain_dx0.1 -mesh ../share/meshes/ex2b_dam_break_dry_domain_dx0.1.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_dry_domain_dx0.1.IC -mannings_n 0 -dt 0.001 -ts_max_time 6 -caseid 3.2
mpiexec -n 4 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_dry_domain_dx0.01 -mesh ../share/meshes/ex2b_dam_break_dry_domain_dx0.01.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_dry_domain_dx0.01.IC -mannings_n 0 -dt 0.001 -ts_max_time 6 -caseid 3.2
mpiexec -n 20 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_dry_domain_dx0.001 -mesh ../share/meshes/ex2b_dam_break_dry_domain_dx0.001.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_dry_domain_dx0.001.IC -mannings_n 0 -dt 0.0001 -ts_max_time 6 -caseid 3.2
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_dam_break_dry_domain_dx0.0001 -mesh ../share/meshes/ex2b_dam_break_dry_domain_dx0.0001.exo -initial_condition ../share/initial_conditions/ex2b_dam_break_dry_domain_dx0.0001.IC -mannings_n 0 -dt 0.00005 -ts_max_time 6 -caseid 3.2
```

Analytical oscillation: parabolic bowl, radially_symmetrical
```
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_radially_symmetrical_T4 -mesh ../share/meshes/ex2b_parabolic_bowl_radially_symmetrical.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_radially_symmetrical.IC -mannings_n 0 -dt 0.001 -ts_max_time 0.56 -caseid 4.1 # T/4
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_radially_symmetrical_T2 -mesh ../share/meshes/ex2b_parabolic_bowl_radially_symmetrical.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_radially_symmetrical.IC -mannings_n 0 -dt 0.001 -ts_max_time 1.121 -caseid 4.1 # T/2
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_radially_symmetrical_T -mesh ../share/meshes/ex2b_parabolic_bowl_radially_symmetrical.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_radially_symmetrical.IC -mannings_n 0 -dt 0.001 -ts_max_time 2.243 -caseid 4.1 # T
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_radially_symmetrical_2T -mesh ../share/meshes/ex2b_parabolic_bowl_radially_symmetrical.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_radially_symmetrical.IC -mannings_n 0 -dt 0.001 -ts_max_time 4.486 -caseid 4.1 # 2T
```

Analytical oscillation: parabolic bowl, planar_surface
```
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_planar_surface_T4 -mesh ../share/meshes/ex2b_parabolic_bowl_planar_surface.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_planar_surface.IC -mannings_n 0 -dt 0.0001 -ts_max_time 1.1214 -caseid 4.2 # T/4
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_planar_surface_T2 -mesh ../share/meshes/ex2b_parabolic_bowl_planar_surface.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_planar_surface.IC -mannings_n 0 -dt 0.0001 -ts_max_time 2.2429 -caseid 4.2 # T/2
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_planar_surface_T -mesh ../share/meshes/ex2b_parabolic_bowl_planar_surface.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_planar_surface.IC -mannings_n 0 -dt 0.0001 -ts_max_time 4.4857 -caseid 4.2 # T
mpiexec -n 40 ./ex2b_validation -savef true -output_prefix ex2b_parabolic_bowl_planar_surface_2T -mesh ../share/meshes/ex2b_parabolic_bowl_planar_surface.exo -initial_condition ../share/initial_conditions/ex2b_parabolic_bowl_planar_surface.IC -mannings_n 0 -dt 0.0001 -ts_max_time 8.9714 -caseid 4.2 # 2T
```

Malpasset Dam break (TELEMAC grid)
```
mpiexec -n 160 ./ex2b_validation -savef true -output_prefix ex2b_malpasset_telemac -mesh ../share/meshes/ex2b_malpasset_telemac.exo -initial_condition ../share/initial_conditions/ex2b_malpasset_telemac.IC -mannings_n 0.033 -dt 0.001 -ts_max_time 4000 -caseid 5.1
```
Malpasset Dam break (PIHM grid)
```
mpiexec -n 160 ./ex2b_validation -savef true -output_prefix ex2b_malpasset_pihm -mesh ../share/meshes/ex2b_malpasset_pihm.exo -initial_condition ../share/initial_conditions/ex2b_malpasset_pihm.IC -mannings_n 0.033 -dt 0.001 -ts_max_time 4000 -caseid 5.2
```
