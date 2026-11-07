# Tutorial 01: ten orifices, 1 m/s

One compression–rebound cycle of the default damper (bore radius 10 mm, rod
radius 4 mm, gas length 120 mm, piston 8 mm, stroke 90 mm, nitrogen at 1 bar
and 300 K) with ten orifices of radius 1 mm at a mean piston speed of 1 m/s,
on the medium mesh, with the k–ω SST model and eight MPI ranks.

```bash
./run.sh check tutorials/01_quick/config.yaml
./run.sh run runs/quick --config tutorials/01_quick/config.yaml
```

`check` prints the derived quantities and the schedule of six mesh segments
on the four sub-meshes `M0, M1, M2, M3, M2, M1`; the stroke reverses inside the
segment on `M3`.

## Expected output

On a workstation with eight cores the run takes about half an hour: about half
a minute for the four sub-meshes (about 107 000 to 91 000 hexahedra), about
25 minutes for the six segments (about 4300 time steps) and a few seconds for
the post-processing and the images. `runs/quick/results/summary.json` then
holds, approximately,

| Key | Value |
|---|---|
| `n_polytropic` | 1.4017 (ratio of specific heats of the gas model: 1.3994) |
| `F_spring_mid_N` | 5.489 |
| `F_damp_mean_N` | 0.031 |
| `E_diss_J` | 0.0057 |
| `Cd_effective` | 0.61 |
| `mass_drift_max_rel` | about 1e-6 |
| `mass_jumps_rel` | five values of order 5e-7 (one per mesh switch) |
| `volume_error_max_rel` | 1.5e-3 (polygonal representation of the curved boundaries) |

The quadrilateral cross-section mesher is not bit-reproducible, so a new run
generates a slightly different cross-section; the spring force and the
polytropic exponent then agree to four digits, the damping force and the
dissipated energy to about 1 %. `results/timeseries.csv` holds one row per
time step, `results/quicklook.png` shows the force loop, and `frames/` the
pressure, velocity and temperature on the plane through the axis and the first
orifice at mid-compression, end of compression and mid-rebound.
