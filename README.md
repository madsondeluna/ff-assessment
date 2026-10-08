# charmmvsgromos

Reproducible comparison of the translation initiation factor eIF4E simulated in
explicit water with an all-atom force field (CHARMM36m) and a united-atom force
field (GROMOS 53A6), using GROMACS 2023 or newer.

The pipeline builds both systems from the same input structure, runs the same
protocol (energy minimization, 100 ps NVT, 100 ps NPT, 100 ns production),
computes the same set of analyses, and produces side-by-side figures and a
summary table.

## Repository layout

```
charmmvsgromos/
├── config.sh                    all parameters (single source of truth)
├── input/                       place the eIF4E structure here (eif4e.pdb)
├── forcefields/                 place the CHARMM36 GROMACS port here
├── mdp/
│   ├── charmm36/                ions, em, nvt, npt, md .mdp files
│   └── gromos53a6/              same files for GROMOS 53A6
├── scripts/
│   ├── lib.sh                   shared bash helpers
│   ├── 00_make_mdps.sh          writes mdp/ from config.sh
│   ├── 01_prepare.sh            build, ions, EM, NVT, NPT, production .tpr
│   ├── 02_run_production.sh     production run (resumes from checkpoint)
│   ├── 03_analysis.sh           RMSD, RMSF, Rg, H-bonds, SASA, mdmat, QC
│   └── 04_plot_comparison.py    comparative figures and statistics
├── requirements.txt             Python dependencies (numpy, matplotlib)
├── runs/                        simulations (created, not tracked)
├── analysis/                    analysis outputs (created, not tracked)
└── figures/                     final figures and summary_statistics.csv
```

## Requirements

| Software | Version | Notes |
|---|---|---|
| GROMACS | 2023 or newer | `gmx` or `gmx_mpi` (set `GMX`); tested with 2023.3 |
| CHARMM36 GROMACS port | February 2026 release recommended | download separately, see below |
| Python | 3.9 or newer | `pip install -r requirements.txt` |
| bash, awk, sed | any recent | standard on Linux and macOS |

GROMOS 53A6 is distributed with GROMACS and needs no download.

### CHARMM36 force field

The CHARMM36 port is not distributed with GROMACS. Download it from the
MacKerell lab (<https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs>) and
extract it into `forcefields/`:

```bash
tar -xzf charmm36-feb2026_cgenff-5.0.ff.tgz -C forcefields/
```

The current ports contain the CHARMM36m protein parameters by default (the
`m` revision of the protein backbone CMAP); lipids, ions and water are the
CHARMM36 ones. If you use another release, set `CHARMM_FF` in `config.sh` to
the directory name without the `.ff` suffix. Do not use the LJ-PME variant of
the port with these `.mdp` files: it requires a different nonbonded protocol.

## Quick start

```bash
# 1. input structure and force field
cp /path/to/eif4e.pdb input/eif4e.pdb
tar -xzf charmm36-feb2026_cgenff-5.0.ff.tgz -C forcefields/

# 2. load GROMACS (example for a module system)
module load gromacs

# 3. build and equilibrate both systems
bash scripts/01_prepare.sh

# 4. production (rerun the same command to resume after an interruption)
bash scripts/02_run_production.sh

# 5. analyses and figures
bash scripts/03_analysis.sh
python3 scripts/04_plot_comparison.py
```

Each script accepts a force field label to process only one system, for
example `bash scripts/02_run_production.sh gromos53a6`. This makes it easy to
run the two production simulations as separate cluster jobs. For queue time
limits use `MAXH=23.5 bash scripts/02_run_production.sh charmm36`; mdrun
stops cleanly and the next call continues from the checkpoint.

Any variable in `config.sh` can be overridden on the command line:

```bash
NT=16 MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu" bash scripts/02_run_production.sh
```

## Configuration

The most relevant variables in `config.sh`:

| Variable | Default | Meaning |
|---|---|---|
| `INPUT_PDB` | `input/eif4e.pdb` | input structure |
| `CHAIN` | `A` | chain(s) to keep; empty string keeps all |
| `HIS_STATES` | `auto` | histidine protonation, one letter per HIS: `D` (ND1), `E` (NE2), `P` (both) |
| `ALLOW_GAPS` | `0` | stop when missing residues are detected |
| `CHARMM_FF` | `charmm36-feb2026_cgenff-5.0` | CHARMM36 port directory in `forcefields/` |
| `BOX_TYPE`, `BOX_DIST` | `dodecahedron`, `1.2` nm | simulation cell |
| `SALT_CONC` | `0.15` M | NaCl added after neutralization |
| `TEMP` | `300` K | temperature |
| `PROD_NS` | `100` | production length (ns) |
| `XTC_PS` | `10` | frame interval (ps) |
| `N_REPLICAS` | `1` | independent replicas per force field |
| `SKIP_NS` | `10` | initial time excluded from averages |
| `NT` | cores minus one | OpenMP threads for mdrun |

The `.mdp` files in `mdp/` are regenerated from `config.sh` every time
`01_prepare.sh` runs. Edit `config.sh` (or `00_make_mdps.sh`), not the
`.mdp` files. The `@SEED@` placeholder in `nvt.mdp` is replaced by a
replica-specific seed (`SEED_BASE + replica`).

## Simulation protocol

The protocol follows the standard lysozyme-in-water workflow. Everything that
defines the protocol is identical in both systems; only the parameters that
belong to each force field's parametrization differ, because applying one
force field's settings to the other is a known source of artefacts and would
bias the comparison.

| Parameter | CHARMM36m | GROMOS 53A6 |
|---|---|---|
| Water model | CHARMM TIP3P (LJ on H) | SPC |
| Hydrogens | all explicit | polar and aromatic only (united-atom CHn) |
| Electrostatics | PME, rc = 1.2 nm | PME, rc = 1.4 nm |
| Lennard-Jones | force-switch 1.0 to 1.2 nm | cutoff 1.4 nm |
| Dispersion correction | none | EnerPres |
| Constraints | bonds to H (LINCS) | all bonds (LINCS) |
| Integrator, dt | leap-frog, 2 fs | leap-frog, 2 fs |
| Energy minimization | steepest descent, Fmax < 1000 kJ/mol/nm | same |
| NVT | 100 ps, V-rescale (0.1 ps), heavy-atom restraints | same |
| NPT | 100 ps, C-rescale (2 ps), heavy-atom restraints | same |
| Production | 100 ns, V-rescale + Parrinello-Rahman, 1 bar | same |
| Coupling groups | Protein, Non-Protein | same |
| Box, ions | dodecahedron 1.2 nm, 0.15 M NaCl | same |

Departures from the original lysozyme tutorial files: C-rescale replaces
Parrinello-Rahman during NPT equilibration (robust far from equilibrium and
available since GROMACS 2021), obsolete options such as `ns_type` are not
used, and genion runs on a `.tpr` built with plain cutoff electrostatics so
that grompp does not stop on the net charge of the unneutralized system.

### Note on GROMOS in modern GROMACS

grompp emits one warning for every GROMOS force field: GROMOS was
parametrized with a twin-range multiple-time-stepping scheme that modern
integrators do not reproduce, so properties such as density can deviate
slightly from the original parametrization (GROMACS issue 2884). The scripts
allow exactly this one extra warning for GROMOS and stop on any other. This
deviation is part of what the comparison measures and should be mentioned in
the methods section.

## Analyses

`03_analysis.sh` first removes periodic boundary artefacts (whole molecules,
protein centred, compact cell), then computes:

| Analysis | GROMACS tool | Atoms | Comparable between AA and UA |
|---|---|---|---|
| RMSD vs minimized structure | `gmx rms` | backbone | yes |
| RMSF per residue | `gmx rmsf -res` | C-alpha | yes |
| Radius of gyration | `gmx gyrate` | protein, mass weighted | yes (united atoms carry the H mass) |
| Intramolecular H-bonds | `gmx hbond` | protein | yes (polar H explicit in both) |
| Protein-water H-bonds | `gmx hbond` | protein, SOL | yes, but includes the water model difference |
| SASA, heavy atoms | `gmx sasa` | protein without H | yes |
| SASA, native | `gmx sasa` | all protein atoms of each model | no, reported for reference only |
| SASA per residue | `gmx sasa -or` | protein without H | yes |
| Contact map | `gmx mdmat` | protein without H | yes |
| Quality checks | `gmx energy` | EM, NVT, NPT, production | n/a |

Hydrogen bonds use the GROMACS geometric criterion (donor-acceptor distance
up to 0.35 nm and hydrogen-donor-acceptor angle up to 30 degrees). GROMACS
2024 rewrote `gmx gyrate` and `gmx hbond`; the script uses the `-legacy`
versions when they exist, so results are produced by the same algorithm on
2023, 2024 and 2025.

Why heavy-atom SASA and contact maps: GROMOS has no aliphatic hydrogens and
`gmx sasa` assigns element radii, so a united-atom CH2 is treated as a bare
carbon and the apolar surface is underestimated. Likewise, explicit hydrogens
shorten residue-residue minimum distances in CHARMM. Removing hydrogens from
both models gives the same atom set with the same radii.

## Outputs

`figures/` (PNG and PDF):

| File | Content |
|---|---|
| `fig1_timeseries` | RMSD, Rg, intramolecular and protein-water H-bonds, SASA |
| `fig2_distributions` | distributions over the equilibrated window |
| `fig3_rmsf` | C-alpha RMSF per residue and CHARMM minus GROMOS difference |
| `fig4_contact_maps` | mean residue minimum distances and their difference |
| `fig5_sasa_per_residue` | heavy-atom SASA per residue and difference |
| `fig6_quality_checks` | EM potential, NVT temperature, NPT pressure and density, production |
| `summary_statistics.csv` | mean, SD and block-averaged SEM per observable |

## Interpreting the comparison

- One 100 ns run per force field shows what happened in that run. To attribute
  a difference to the force field rather than to run-to-run variability, use
  `N_REPLICAS=3` or more; per-residue figures then show replica means with
  SD bands, and the table reports the SD of the replica means.
- Consecutive frames are correlated. The table reports block-averaged SEM
  (5 blocks) for single replicas; frame-level tests (for example a t-test over
  all frames) overstate significance.
- The water model is part of each force field (SPC for GROMOS, TIP3P for
  CHARMM). Differences in protein-water H-bonds reflect both the protein and
  the water parameters.
- Histidine protonation and net charge are cross-checked after `pdb2gmx`; the
  script warns if they differ between force fields. Set `HIS_STATES` to fix
  them explicitly.
- Crystal structures of eIF4E usually lack the flexible N-terminal region.
  Internal gaps stop the preparation (model them first); a truncated terminus
  is accepted and receives charged termini in both force fields.
- Disk use is on the order of 1 to 2 GB per 100 ns replica, dominated by the
  full-system trajectory needed for protein-water H-bonds.

## References

- Huang J. et al. CHARMM36m: an improved force field for folded and
  intrinsically disordered proteins. Nat. Methods 14, 71-73 (2017).
- Oostenbrink C. et al. A biomolecular force field based on the free enthalpy
  of hydration and solvation: the GROMOS force-field parameter sets 53A5 and
  53A6. J. Comput. Chem. 25, 1656-1676 (2004).
- Abraham M. J. et al. GROMACS: high performance molecular simulations
  through multi-level parallelism from laptops to supercomputers. SoftwareX
  1-2, 19-25 (2015).
- Lemkul J. A. From proteins to perturbed Hamiltonians: a suite of tutorials
  for the GROMACS-2018 molecular simulation package. Living J. Comput. Mol.
  Sci. 1, 5068 (2018).

## Author

Madson A. de Luna Aragão
