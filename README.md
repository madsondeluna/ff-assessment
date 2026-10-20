# ff-assessment

Reproducible four-way comparison of a small protein simulated in explicit
water with two all-atom force fields (CHARMM36m and AMBER ff14SB) and two
united-atom force fields (GROMOS 53A6 and GROMOS 54A7), using GROMACS 2023 or
newer. Each force field is run in three independent replicas (12 production
simulations of 100 ns by default), strictly one after the other, with a
progress log.

The default test protein is human ubiquitin (PDB 1UBQ, 76 residues, 1.8 A):
one alpha helix, one 3-10 helix and a five-stranded mixed beta sheet, no
gaps and no disulfides, which makes it a good small alpha/beta case for
asking how each force field keeps (or loses) secondary structure. The
structure is downloaded from the RCSB and cleaned automatically; any other
entry can be used with `PDB_ID`, or any prepared structure with `INPUT_PDB`.

The pipeline builds every system from the same input structure, runs the same
protocol (energy minimization, 100 ps NVT, 100 ps NPT, 100 ns production),
computes the same set of analyses for every replica, and produces two levels
of comparison, as in a multi-replica benchmark:

- between replicas of the same force field (convergence, run-to-run spread);
- between force fields (replica means with SD, pooled distributions, and
  "force field minus reference" panels, with CHARMM36m as the default
  reference), including secondary structure per residue and time (DSSP:
  helix dark blue, strand red, coil white).

The repository started as a two-way CHARMM36m vs GROMOS 53A6 comparison
(`charmmvsgromos`) and was extended to GROMOS 54A7 and AMBER; the old name
still redirects to this repository on GitHub.

## Repository layout

```
ff-assessment/
├── config.sh                    all parameters (single source of truth)
├── input/                       downloaded and cleaned structure (1ubq.pdb)
├── forcefields/                 place downloaded force field ports here
├── mdp/
│   ├── charmm36m/               ions, em, nvt, npt, md .mdp files
│   ├── amber/                   same files for AMBER
│   ├── gromos53a6/              same files for GROMOS 53A6
│   └── gromos54a7/              same files for GROMOS 54A7
├── scripts/
│   ├── lib.sh                   shared bash helpers (force field table, PDB cleaning)
│   ├── check.sh                 pre-flight check: GROMACS, force fields, Python, disk
│   ├── get_structure.sh         download PDB_ID from the RCSB and clean it
│   ├── 00_make_mdps.sh          writes mdp/ from config.sh
│   ├── 01_prepare.sh            build, ions, EM, NVT, NPT, production .tpr
│   ├── 02_run_production.sh     production runs, one at a time (resume from checkpoint)
│   ├── 03_analysis.sh           RMSD, RMSF, Rg, H-bonds, SASA, mdmat, DSSP, QC
│   ├── 04_plot_comparison.py    comparative figures and statistics
│   ├── run_all.sh               the whole campaign in sequence, logged
│   └── status.sh                progress table of every run
├── requirements.txt             Python dependencies (numpy, matplotlib; scipy optional)
├── logs/                        pipeline and production logs (created, not tracked)
├── runs/                        simulations (created, not tracked)
├── analysis/                    analysis outputs (created, not tracked)
└── figures/                     final figures, summary_statistics.csv, summary_table.csv
```

## Setup on a new machine (Ubuntu)

Five commands, then a check that tells you what is still missing and how to
get it:

```bash
sudo apt install gromacs python3-venv curl          # GROMACS 2023.3 on Ubuntu 24.04
git clone https://github.com/madsondeluna/ff-assessment.git
cd ff-assessment
python3 -m venv ~/.venvs/ffa && source ~/.venvs/ffa/bin/activate && pip install -r requirements.txt
bash scripts/check.sh
```

`check.sh` finds the GROMACS executable on its own (the first of `gmx`,
`gmx_mpi`, `gmx_d`, `gmx_mpi_d` in `PATH`; on a machine with a single
installation, from apt, conda or a module, nothing needs to be configured),
prints its version and whether it has GPU support, and then lists every force
field with `[ok]` or `[MISSING]`. For a missing one it prints the download
address and the exact commands; GROMOS 53A6, GROMOS 54A7 and the AMBER
fallback ship with GROMACS, so on a fresh machine the only thing it asks for
is the CHARMM36 port (next section). `run_all.sh` runs the same check first
and stops if something required is absent. Example of the output on a
machine without the port:

```
Force fields (FF_LABELS: charmm36m amber gromos53a6 gromos54a7)
  [MISSING] charmm36m: directory 'charmm36-feb2026_cgenff-5.0.ff' not found in forcefields/ ...
            The CHARMM36 port is not distributed with GROMACS. Download it by hand:
              1. open https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs in a browser
              ...
  [ok]      amber: AMBER ff99SB-ILDN -> /usr/share/gromacs/top/amber99sb-ildn.ff (water: tip3p)
  [note]    ff14SB is not in this GROMACS (built in from 2026); using ff99SB-ILDN ...
  [ok]      gromos53a6: GROMOS 53A6 -> /usr/share/gromacs/top/gromos53a6.ff (water: spc)
  [ok]      gromos54a7: GROMOS 54A7 -> /usr/share/gromacs/top/gromos54a7.ff (water: spc)
```

## Requirements

| Software | Version | Notes |
|---|---|---|
| GROMACS | 2023 or newer | autodetected; set `GMX` only when several versions are installed; 2026 or newer for the built-in AMBER ff14SB/ff19SB |
| CHARMM36 GROMACS port | most recent release | downloaded by hand from the MacKerell lab, see below |
| Python | 3.9 or newer | `pip install -r requirements.txt` (numpy, matplotlib; scipy optional) |
| bash, awk, sed, curl or wget | any recent | standard on Linux and macOS |

GROMOS 53A6, GROMOS 54A7 and AMBER ff99SB-ILDN are distributed with every
GROMACS version and need no download. A GPU is optional: without one the
runs are slower but identical.

### Force field labels

Every force field is identified by a label used for directory names and
plots. The label maps to a force field directory, a water model and a family;
the family decides the nonbonded and constraint settings.

| Label | Family | Force field directory | Water | Where it comes from |
|---|---|---|---|---|
| `charmm36m` | charmm | `CHARMM_FF` (default `charmm36-feb2026_cgenff-5.0`) | CHARMM TIP3P | MacKerell lab port, downloaded by hand into `forcefields/` |
| `amber` | amber | `AMBER_FF` (default `auto`: ff14SB if available, else ff99SB-ILDN) | `AMBER_WATER` (default TIP3P) | built into GROMACS, or a port in `forcefields/` |
| `gromos53a6` | gromos | `gromos53a6` | SPC | built into GROMACS |
| `gromos54a7` | gromos | `gromos54a7` | SPC | built into GROMACS |

### CHARMM36m: manual download

The CHARMM36 port is not distributed with GROMACS and cannot be fetched by a
script (the site asks for a registration). Do it once by hand:

1. Open <https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs> in a
   browser and register with an e-mail address.
2. Under "GROMACS", download the most recent `charmm36-<month><year>.ff.tgz`
   (for example `charmm36-feb2026_cgenff-5.0.ff.tgz`). Do not take the
   `lj-pme` variant: it needs a different nonbonded protocol.
3. Copy the file to the machine that will run the simulations and extract it
   into `forcefields/`:

   ```bash
   scp ~/Downloads/charmm36-*.ff.tgz user@workstation:~/ff-assessment/forcefields/
   tar -xzf forcefields/charmm36-*.ff.tgz -C forcefields/
   ls forcefields/            # one directory ending in .ff
   ```

4. If the directory is not `charmm36-feb2026_cgenff-5.0.ff`, tell the scripts
   its name without the `.ff` suffix, either on the command line or by
   editing `CHARMM_FF` in `config.sh`:

   ```bash
   export CHARMM_FF=charmm36-jul2022
   ```

5. Run `bash scripts/check.sh` again: the `charmm36m` line must now say
   `[ok]`.

The current ports contain the CHARMM36m protein parameters by default (the
`m` revision of the protein backbone CMAP); lipids, ions and water are the
CHARMM36 ones. To start the other three force fields while the port is on
its way: `FF_LABELS="amber gromos53a6 gromos54a7" bash scripts/run_all.sh`;
running the script again later with the default labels adds the CHARMM runs
without touching the finished ones.

### AMBER

The label `amber` is one AMBER protein force field, chosen with `AMBER_FF`
and `AMBER_WATER`:

| `AMBER_FF` | `AMBER_WATER` | Availability | Notes |
|---|---|---|---|
| `auto` (default) | `tip3p` | any | `amber14sb` when the installed GROMACS (or `forcefields/`) has it, otherwise `amber99sb-ildn`, with a warning in the log and the name in the figures |
| `amber14sb` | `tip3p` | GROMACS >= 2026 | ff14SB, the usual reference AMBER force field; 3-site water like the other three systems |
| `amber19sb` | `opc` (or `opc3`) | GROMACS >= 2026 | ff19SB with per-residue CMAP; OPC is a 4-site model (about 30 percent more water sites). The release notes call both ports validation pending |
| `amber99sb-ildn` | `tip3p` | every GROMACS | what `auto` picks on GROMACS 2023 to 2025 (Ubuntu 24.04 ships 2023.3) |

Older GROMACS versions can also use a user-contributed ff14SB/ff19SB port
extracted into `forcefields/` and named in `AMBER_FF`, as for CHARMM (the
`check.sh` output says so when ff14SB is not found). The
scripts check that ff19SB is paired with OPC and that the solvent box used
by `gmx solvate` matches the water model (3-site or 4-site). The ion
molecule names are read from the ion `.itp` that `pdb2gmx` includes in the
topology, which in the built-in ff14SB/ff19SB ports depends on the water
model.

## Quick start

```bash
# 0. (module systems only) load GROMACS; apt and conda installs are found by themselves
module load gromacs

# 1. check: GROMACS found, force fields found, Python packages installed
bash scripts/check.sh

# 2. optional: a 10 minute smoke test of the whole pipeline before the real campaign
PROD_NS=0.1 N_REPLICAS=1 NVT_PS=10 NPT_PS=10 bash scripts/run_all.sh && rm -rf runs analysis figures logs

# 3. everything, in sequence, detached, with a log (GPU: see below)
nohup bash scripts/run_all.sh > /dev/null 2>&1 &

# 4. follow the campaign
tail -f logs/pipeline_latest.log      # every line of every stage
bash scripts/status.sh                # one line per run: stage, ns done, %, ns/day
```

With an NVIDIA GPU and a CUDA build of GROMACS (`check.sh` prints "GPU
support: CUDA"), export `MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu"` before
step 3. With the Ubuntu apt build (CPU only) leave it empty.

`run_all.sh` chains `check.sh`, `get_structure.sh` (download and clean
1UBQ), `01_prepare.sh` (four systems, three replicas each),
`02_run_production.sh` (12 runs), `03_analysis.sh` and
`04_plot_comparison.py`; each stage starts only when the previous one has
finished, and rerunning the script skips what is already done. The stages
can also be run one by one:

```bash
bash scripts/get_structure.sh          # input/1ubq.pdb (and 1ubq_raw.pdb)
bash scripts/01_prepare.sh
nohup bash scripts/02_run_production.sh > logs/production_stdout.log 2>&1 &
bash scripts/03_analysis.sh
python3 scripts/04_plot_comparison.py
```

Each script accepts force field labels to process a subset, for example
`bash scripts/02_run_production.sh gromos54a7 amber`, and `REPS="2 3"`
selects replicas, so the production runs can be split into separate jobs. For
queue time limits use `MAXH=23.5 bash scripts/02_run_production.sh charmm36m`;
mdrun stops cleanly and the next call continues from the checkpoint.

### Sequential execution and logs

Production runs are strictly sequential: `02_run_production.sh` starts the
next `mdrun` only when the previous one has finished (or stopped at `MAXH`),
so one simulation has all threads and the GPU to itself. `RUN_ORDER=replica`
(default) runs replica 1 of every force field first, then replica 2 and 3,
so that a one-replica comparison is available after the first four runs;
`RUN_ORDER=forcefield` finishes one force field before starting the next.

| Log | Content |
|---|---|
| `logs/pipeline_<timestamp>.log` (`pipeline_latest.log` links to the newest) | everything printed by `run_all.sh` and the stages it calls |
| `logs/production.log` | one line per run start, resume and finish, with wall time and ns/day |
| `runs/<ff>/rep<N>/md_stdout.log`, `md.log` | mdrun output of each run |
| `bash scripts/status.sh` | table: build and equilibration state, production state (running, stopped, finished), ns done, percent, ns/day, analysis up to date or stale |

Any variable in `config.sh` can be overridden on the command line or
exported in the shell before `nohup`:

```bash
NT=9 MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu" bash scripts/02_run_production.sh
GMX=/opt/gromacs-2026/bin/gmx bash scripts/check.sh            # several GROMACS versions installed
AMBER_FF=amber19sb AMBER_WATER=opc bash scripts/01_prepare.sh amber
FF_LABELS="amber gromos53a6 gromos54a7" bash scripts/run_all.sh  # without CHARMM for now
PDB_ID=1PGA INPUT_PDB=input/1pga.pdb bash scripts/run_all.sh   # another protein
PROD_NS=5 N_REPLICAS=1 bash scripts/run_all.sh                # short test
```

Ubiquitin in a 1.2 nm dodecahedron is a small system (around 20 000 atoms
all-atom, fewer united-atom), so with one GPU each 100 ns run takes hours
rather than days and the 12 runs complete in a few days; the united-atom
systems are the fastest.

## Configuration

The most relevant variables in `config.sh`:

| Variable | Default | Meaning |
|---|---|---|
| `GMX` | empty (autodetect) | GROMACS executable; set only when several versions are installed |
| `PDB_ID` | `1UBQ` | entry downloaded and cleaned when `INPUT_PDB` is missing (empty disables) |
| `INPUT_PDB` | `input/1ubq.pdb` | input structure |
| `CHAIN` | `A` | chain(s) to keep; empty string keeps all |
| `HIS_STATES` | `auto` | histidine protonation, one letter per HIS: `D` (ND1), `E` (NE2), `P` (both) |
| `ALLOW_GAPS` | `0` | stop when missing residues are detected |
| `FF_LABELS` | `charmm36m amber gromos53a6 gromos54a7` | force fields to build and compare (a space-separated string when set from the environment) |
| `FF_REFERENCE` | `charmm36m` | reference for the difference panels (also `--reference` in the plot script) |
| `CHARMM_FF` | `charmm36-feb2026_cgenff-5.0` | CHARMM36 port directory in `forcefields/` |
| `AMBER_FF`, `AMBER_WATER` | `auto`, `tip3p` | AMBER force field (ff14SB if available, else ff99SB-ILDN) and water model |
| `BOX_TYPE`, `BOX_DIST` | `dodecahedron`, `1.2` nm | simulation cell |
| `SALT_CONC` | `0.15` M | NaCl added after neutralization |
| `TEMP` | `300` K | temperature |
| `PROD_NS` | `100` | production length (ns) |
| `XTC_PS` | `10` | frame interval (ps) |
| `N_REPLICAS` | `3` | independent replicas per force field |
| `RUN_ORDER` | `replica` | sequence of the production runs (`replica` or `forcefield`) |
| `SKIP_NS` | `10` | initial time excluded from averages |
| `DSSP_DT_PS` | `100` | frame interval for the DSSP assignment (ps) |
| `NT` | cores minus one | OpenMP threads for mdrun |

The `.mdp` files in `mdp/` are regenerated from `config.sh` every time
`01_prepare.sh` runs. Edit `config.sh` (or `00_make_mdps.sh`), not the
`.mdp` files. The `@SEED@` placeholder in `nvt.mdp` is replaced by a
replica-specific seed (`SEED_BASE + replica`), so every replica starts the
NVT equilibration from the minimized structure with its own velocities and
goes through its own NPT equilibration before production.

## Simulation protocol

The protocol follows the standard lysozyme-in-water workflow. Everything that
defines the protocol is identical in every system; only the parameters that
belong to each force field family's parametrization differ, because applying
one family's settings to another is a known source of artefacts and would
bias the comparison.

| Parameter | CHARMM36m | AMBER ff14SB | GROMOS 53A6 | GROMOS 54A7 |
|---|---|---|---|---|
| Model | all-atom | all-atom | united-atom (CHn) | united-atom (CHn) |
| Water model | CHARMM TIP3P (LJ on H) | TIP3P | SPC | SPC |
| Hydrogens | all explicit | all explicit | polar and aromatic only | polar and aromatic only |
| Electrostatics | PME, rc = 1.2 nm | PME, rc = 1.0 nm | PME, rc = 1.4 nm | PME, rc = 1.4 nm |
| Lennard-Jones | force-switch 1.0 to 1.2 nm | cutoff 1.0 nm, potential shift | cutoff 1.4 nm, potential shift | cutoff 1.4 nm, potential shift |
| Dispersion correction | none | EnerPres | EnerPres | EnerPres |
| Constraints | bonds to H (LINCS) | bonds to H (LINCS) | all bonds (LINCS) | all bonds (LINCS) |
| Integrator, dt | leap-frog, 2 fs | same | same | same |
| Energy minimization | steepest descent, Fmax < 1000 kJ/mol/nm | same | same | same |
| NVT | 100 ps, V-rescale (0.1 ps), heavy-atom restraints | same | same | same |
| NPT | 100 ps, C-rescale (2 ps), heavy-atom restraints | same | same | same |
| Production | 100 ns, V-rescale + Parrinello-Rahman, 1 bar | same | same | same |
| Coupling groups | Protein, Non-Protein | same | same | same |
| Box, ions | dodecahedron 1.2 nm, 0.15 M NaCl | same | same | same |
| Replicas | 3 (seeds 2027, 2028, 2029) | same | same | same |

GROMOS 53A6 and 54A7 share the same functional form and nonbonded protocol;
54A7 is a reparametrization (backbone dihedrals, charged side chains, an
extra atom type) that corrects the helix propensity of 53A6, so running both
isolates the effect of the parameters from the effect of the united-atom
model. AMBER ff14SB uses the rc = 1.0 nm cutoff with a long-range dispersion
correction, which is the usual GROMACS setup for AMBER force fields (the
family was parametrized with an 8 to 10 Angstrom cutoff plus PME). With
`AMBER_FF=amber19sb` the same block applies; the water model changes to OPC.

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
allow exactly this one extra warning for the two GROMOS labels and stop on
any other. This deviation is part of what the comparison measures and should
be mentioned in the methods section.

## Analyses

`03_analysis.sh` first removes periodic boundary artefacts (whole molecules,
protein centred, compact cell), then computes, for every replica:

| Analysis | GROMACS tool | Atoms | Comparable between AA and UA |
|---|---|---|---|
| RMSD vs minimized structure | `gmx rms` | backbone | yes |
| RMSF per residue | `gmx rmsf -res` | C-alpha | yes |
| Radius of gyration | `gmx gyrate` | protein, mass weighted | yes (united atoms carry the H mass) |
| Intramolecular H-bonds | `gmx hbond` | protein | yes (polar H explicit in all) |
| Protein-water H-bonds | `gmx hbond` | protein, SOL | yes, but includes the water model difference |
| SASA, heavy atoms | `gmx sasa` | protein without H | yes |
| SASA, native | `gmx sasa` | all protein atoms of each model | no, reported for reference only |
| SASA per residue | `gmx sasa -or` | protein without H | yes |
| Contact map | `gmx mdmat` | protein without H | yes |
| Secondary structure | `gmx dssp -hmode dssp` | backbone N, CA, C, O | yes (amide H rebuilt from the backbone geometry) |
| Quality checks | `gmx energy` | EM, NVT, NPT, production | n/a |

Hydrogen bonds use the GROMACS geometric criterion (donor-acceptor distance
up to 0.35 nm and hydrogen-donor-acceptor angle up to 30 degrees). GROMACS
2024 rewrote `gmx gyrate` and `gmx hbond`; the script uses the `-legacy`
versions when they exist, so results are produced by the same algorithm on
2023 and newer.

Why heavy-atom SASA and contact maps: GROMOS has no aliphatic hydrogens and
`gmx sasa` assigns element radii, so a united-atom CH2 is treated as a bare
carbon and the apolar surface is underestimated. Likewise, explicit hydrogens
shorten residue-residue minimum distances in CHARMM and AMBER. Removing
hydrogens from every model gives the same atom set with the same radii.

Secondary structure: `gmx dssp` (native since GROMACS 2023) assigns the
eight DSSP states per residue and frame. `-hmode dssp` rebuilds the amide
hydrogen from the C and O of the previous residue, as the original DSSP
does, so the assignment uses the same backbone geometry rule for every force
field instead of each model's own hydrogen positions. The plot script reduces
the states to three classes, as usual in force field benchmarks: helix (H,
G, I), strand (E, B) and coil (T, S, P, loop). The input structure is
assigned with the same tool and used as the reference for the per-residue
panels and for `ss_retained`.

## Outputs

`figures/` (PNG and PDF):

| File | Content |
|---|---|
| `fig1_timeseries` | RMSD, Rg, intramolecular and protein-water H-bonds, SASA: replica mean with min-max band per force field |
| `fig2_distributions` | distributions over the equilibrated window, pooled over replicas |
| `fig3_rmsf` | C-alpha RMSF per residue (replica mean, SD band) and difference of each force field vs the reference |
| `fig4_contact_maps` | mean residue minimum distances per force field and differences vs the reference |
| `fig5_sasa_per_residue` | heavy-atom SASA per residue and differences vs the reference |
| `fig6_quality_checks` | EM potential, NVT temperature, NPT pressure and density, production (line style: replica) |
| `fig7_replicas` | one row per force field, each replica as its own line: RMSD, Rg, H-bonds, SASA, RMSF |
| `fig8_dssp_maps` | secondary structure per residue (y) and time (x), one panel per replica, rows per force field; helix dark blue, strand red, coil white |
| `fig9_dssp_content` | helix and strand content vs time (replica mean, min-max band, input structure as dashed line) and per-residue helix and strand probability over the averaging window against the input structure (hatched) |
| `summary_statistics.csv` | mean, SD and block-averaged SEM per observable, replica and force field |
| `summary_table.csv` | one row per force field and observable: mean, SD over replicas, difference vs reference, Welch's t-test |

`summary_table.csv` is the table to report (mean plus or minus SD of the
three replica means for RMSD, Rg, mean RMSF, H-bonds, SASA, percent of
residues in helix, strand and coil, and `ss_retained`, the percent of
residues that keep the class they have in the input structure), and the
console prints it at the end of the plot script. `--reference` changes the
reference force field and `--only` restricts the figures to a subset of
labels.

## Interpreting the comparison

- Replica spread first, force field differences second. `fig7_replicas` and
  the SD column show whether three runs of the same force field agree; a
  difference between force fields that is smaller than the replica spread is
  not a force field effect. With n = 3 the Welch test has 2 to 4 degrees of
  freedom, so only large and consistent differences reach p < 0.05; report
  the difference and its uncertainty rather than the p-value alone.
- Consecutive frames are correlated. `summary_statistics.csv` reports a
  block-averaged SEM (5 blocks) per replica; frame-level tests (for example
  a t-test over all frames) overstate significance.
- The water model is part of each force field (SPC for GROMOS, CHARMM TIP3P
  for CHARMM36m, TIP3P or OPC for AMBER). Differences in protein-water
  H-bonds reflect both the protein and the water parameters.
- GROMOS 53A6 vs 54A7 isolates a reparametrization within one family;
  CHARMM36m vs AMBER compares two all-atom families; all-atom vs GROMOS
  compares resolutions. Keep these three contrasts separate when reading the
  difference panels.
- Secondary structure is where the families are expected to differ most:
  53A6 is known to under-stabilize helices (the reason for 54A7), and the
  all-atom families differ in helix propensity. `fig8_dssp_maps` shows
  whether a loss is a transient fraying of the termini of a helix or strand
  or a permanent unfolding, and whether it happens in every replica;
  `ss_retained` summarizes it per run.
- Histidine protonation and net charge are cross-checked after `pdb2gmx` for
  every force field against the first one built; the script warns if they
  differ. Set `HIS_STATES` to fix them explicitly.
- 1UBQ is complete (residues 1 to 76); its C-terminal tail (72 to 76) is
  flexible and dominates the RMSF. For another entry, internal gaps stop the
  preparation (model them first); a truncated terminus is accepted and
  receives charged termini in every force field.
- Disk use is below 1 GB per 100 ns replica for ubiquitin (around 20 000
  atoms, a frame every 10 ps), dominated by the full-system trajectory
  needed for protein-water H-bonds.

## Migrating from the two-force-field version

The label `charmm36` was renamed `charmm36m`. Existing runs can be kept with
`mv runs/charmm36 runs/charmm36m && mv analysis/charmm36 analysis/charmm36m`;
the default number of replicas changed from 1 to 3, so `01_prepare.sh` will
add `rep2` and `rep3` next to an existing `rep1` without touching it. The
default protein changed from eIF4E to ubiquitin; to keep eIF4E, set
`PDB_ID=""` and `INPUT_PDB=input/eif4e.pdb`. Existing analyses lack
`dssp.dat`; `FORCE=1 bash scripts/03_analysis.sh` recomputes them.

## References

- Huang J. et al. CHARMM36m: an improved force field for folded and
  intrinsically disordered proteins. Nat. Methods 14, 71-73 (2017).
- Maier J. A. et al. ff14SB: improving the accuracy of protein side chain and
  backbone parameters from ff99SB. J. Chem. Theory Comput. 11, 3696-3713
  (2015).
- Tian C. et al. ff19SB: amino-acid-specific protein backbone parameters
  trained against quantum mechanics energy surfaces in solution. J. Chem.
  Theory Comput. 16, 528-552 (2020).
- Lindorff-Larsen K. et al. Improved side-chain torsion potentials for the
  Amber ff99SB protein force field. Proteins 78, 1950-1958 (2010).
- Oostenbrink C. et al. A biomolecular force field based on the free enthalpy
  of hydration and solvation: the GROMOS force-field parameter sets 53A5 and
  53A6. J. Comput. Chem. 25, 1656-1676 (2004).
- Schmid N. et al. Definition and testing of the GROMOS force-field versions
  54A7 and 54B7. Eur. Biophys. J. 40, 843-856 (2011).
- Abraham M. J. et al. GROMACS: high performance molecular simulations
  through multi-level parallelism from laptops to supercomputers. SoftwareX
  1-2, 19-25 (2015).
- Lemkul J. A. From proteins to perturbed Hamiltonians: a suite of tutorials
  for the GROMACS-2018 molecular simulation package. Living J. Comput. Mol.
  Sci. 1, 5068 (2018).
- Kabsch W., Sander C. Dictionary of protein secondary structure: pattern
  recognition of hydrogen-bonded and geometrical features. Biopolymers 22,
  2577-2637 (1983).
- Vijay-Kumar S., Bugg C. E., Cook W. J. Structure of ubiquitin refined at
  1.8 A resolution. J. Mol. Biol. 194, 531-544 (1987).

## Author

Madson A. de Luna Aragão
