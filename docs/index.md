# SlimFlow Documentation

SlimFlow is a small Python module that makes running and processing [SLiM
simulations](https://messerlab.org/slim/) with
[Snakemake](https://snakemake.readthedocs.io/en/stable/) much easier. It will
also work with [msprime](https://tskit.dev/software/msprime.html), though this
is currently in development. I've used this module for several projects (e.g.
Buffalo and Coop 2020, Buffalo and Kern 2024), and split it out so others can
use. 

SlimFlow takes a simple YAML configuration file of simulation parameters and
generates all possible parameter combinations with their replicates, and
uniquely seeds each one, creating a list of *filepath targets*. SlimFlow also
handles (1) creating the *filepath templates*, which maps simulation parameters
to Snakemake wildcards, and (2) passing parameters to SLiM through the command
line. The end result is that fairly complicated simulations can be created and
executed via Snakemake with very little boilerplate code. Below is an example
of running some BGS simulations according to the arbitrary parameters specified
in a YAML file:

```python
from slimflow import GridRuns

if not len(config):
  raise ValueError("config file not specified, use --configfile <config>.yml")

run = GridRuns(config)

# generate all the targets for this YAML file.
all_sims = run.generate_targets()

# a simple, small amount of boilerplate code will run all simulations.
rule slim:
  input: run.script, **run.input
  output: **run.output_template()
  shell: run.slim_cmd()

rule all:
  input: all_sims['filepath']

```

The `GridRuns` class is designed to generate the output files in a
well-organized fashion, in a tidy readable directory structure that is easy 
to process like:

```
<name>/<dir>/<key1__val1>/<key2__val1>/rep_<rep>__seed_<seed>__<suffix>
```

Here is a concrete example of a file structure from the BGS simulation example
in `examples/bgs`.

```
$ tree runs
runs
└── bgs
    └── N__1000
        ├── mu__1e-08
        │   ├── sh__0.01
        │   │   └── rbp__1e-08
        │   │       ├── rep_0__seed_7919168045412322065__log.tsv.gz
        │   │       ├── rep_0__seed_7919168045412322065__summary.tsv
        │   │       ├── rep_0__seed_7919168045412322065__treeseq.tree
        │   │       ├── rep_1__seed_7919168045412322066__log.tsv.gz
        │   │       ├── rep_1__seed_7919168045412322066__summary.tsv
        │   │       └── rep_1__seed_7919168045412322066__treeseq.tree
        │   └── sh__0.1
        │       └── rbp__1e-08
        │           ├── rep_0__seed_6432084778622665797__log.tsv.gz
        │           ├── rep_0__seed_6432084778622665797__summary.tsv
[...]
```

SlimFlow generates all combinations of the variables, and creates sequential
random seeds per parameter combination for each seed. This way, if more
replicates are added, they do alter the state of the pseudorandom number
generator, which would trigger a bunch of unnecessary simulations.

SlimFlow is designed to simplify writing Snakemake rules for running such SLiM
simulations. It will automatically generate list of *target files* (e.g. the
resulting files of a simulation) for all parameter combinations and replicates,
as well as *target templates* which wildcard templates that match the file
pattern template of these target files. This may seem abstract, but it
eliminates a great deal of boilerplate Snakemake stuff.

## Example

The best way to see how SlimFlow works is through a simple example; the full
source of the example is in `examples/bgs`. Suppose you want to run some
background selection simulation with a SLiM file `bgs.slim`, across several
different combinations of parameters and independent replicates with different
seeds. 

### YAML Simulation Parameter Configuration

First, we would write a YAML config file for this simulation "run" that looks
like:

```yaml
---
name: bgs
script: bgs.slim
dir: runs
nreps: 2
seed: 42
suffices:
  treeseq_file: treeseq.tree
  log_file: log.tsv.gz
variables:
  N:
  - 1000
  mu:
  - 2.0e-08
  - 1.0e-08
  sh:
  - 0.01
  - 0.1
  rbp:
  - 1.0e-08
```

This specifies the key bits of a bunch of simulation runs: the master `seed`,
the number of replicates (`nreps`), the simulation results directory (`dir`),
and most importantly, the parameters *variables* that are passed to the SLiM
simulation as command line arguments.

In some cases, simulations have multiple output files, such as a tree sequence
`.tree` file, as well as a `.tsv.gz` log file. SlimFlow supports that through
different *suffices* that will have the same output filename, but different
suffices. This common filename up to the suffix allows them to be paired in
downstream analyses. 

### Generating Simulation Targets

Then, using SlimFlow, we can load in the config YAML file:

```python
>>> import yaml
>>> from slimflow import GridRuns
>>> config = yaml.safe_load(open('examples/bgs/bgs_config.yml'))
>>> run = GridRuns(config)
>>> run
GridRuns(bgs):
Number of parameter combinations (𝑘): 4
Number of replicates (𝑛): 2
Total simulations (𝑘 × 𝑛): 8
Variables:
  N ∈ {1000}
  mu ∈ {2e-08, 1e-08}
  sh ∈ {0.01, 0.1}
  rbp ∈ {1e-08}
```

All simulation run output files (i.e. "targets") can be generated with
`GridRuns.generate_targets()`. This method just generates all parameter
combinations and seeds them, creating the *seeded filename* which looks like
`rep_0__seed_7138484576005690179__<suffix>`. The suffices for simulation output
files (e.g. `treeseq.tree` or `log.tsv.gz`) are defined in the YAML
configuration file; these are the files expected to result from each SLiM
simulation.

All the results are stored in [Polars](https://pola.rs) dataframe, which
contain a `filepath` column to each simulation file, as well as each parameter
combination. Here's the first few rows and a selection of columns:

```python
>>> run.generate_targets().select(['N', 'mu', 'sh', 'filepath']).head(2)
shape: (2, 4)
┌──────┬───────────┬──────┬───────────────────────────────────┐
│ N    ┆ mu        ┆ sh   ┆ filepath                          │
│ ---  ┆ ---       ┆ ---  ┆ ---                               │
│ i64  ┆ f64       ┆ f64  ┆ str                               │
╞══════╪═══════════╪══════╪═══════════════════════════════════╡
│ 1000 ┆ 2.0000e-8 ┆ 0.01 ┆ runs/bgs/N__1000/mu__2e-08/sh__0… │
│ 1000 ┆ 2.0000e-8 ┆ 0.01 ┆ runs/bgs/N__1000/mu__2e-08/sh__0… │
└──────┴───────────┴──────┴───────────────────────────────────┘
```

This dataframe can then be used to join in results during downstream
processing. Here is an example of generating target files with the default
suffices:

```python
>>> targets = run.generate_targets()
# targets were generated for both suffices for this simulation replicate:
>>> targets[0]
'runs/region/N__1000/mu__2e-08/sh__0.01/rbp__1e-08/rep0_seed7138484576005690179_treeseq.tree'
>>> targets[1]
'runs/region/N__1000/mu__2e-08/sh__0.01/rbp__1e-08/rep0_seed7138484576005690179_log.tsv.gz'
```

### Generating Snakefile Input/Output Templates

Additionally, SlimFlow generates a dictionary containing the filename tempate, or
pattern with wildcards that are matched by rules so parameters can flow into
the `run` or `shell` directives.

These templates can be directly used in a Snakefile rule's `input` and `output`
blocks:

```python
>>> print(run.output_template())
{'treeseq_file': 'runs/region/N__{N}/mu__{mu}/sh__{sh}/rbp__{rbp}/rep{rep}_seed{seed}_treeseq.tree', 
 'log_file': 'runs/region/N__{N}/mu__{mu}/sh__{sh}/rbp__{rbp}/rep{rep}_seed{seed}_log.tsv.gz'} 
```

Note that suffices are stored in dictionaries so that multiple input/output
files could be used with rules. These dictionaries are unpacked with the `**`
before run into `key=value` arguments for `input`/`output`.

Both `GridRuns.generate_targets()` and `GridRuns.output_template()` can also create
target files with different suffices and parent output directories (i.e. like
`runs` in the example config). This is to make downstream analyses easier: for
example, if one Snakemake rule needed to take in tree sequences and add
mutations, it could use `input: **run.output_template(suffices={'tree':
'treeseq.tree'})` and `output: **run.output_template(suffices={'tree':
'mutations.tree'})` to create these target templates. Similarly, all of the
target files could be created with `run.generate_targets(suffices={'tree':
'mutations.tree'})`. 

### Generation of SLiM Commands for Snakemake Rules

SlimFlow also generates the SLiM command string used to pass in the parameter
arguments as wildcards:

```python
>>> print(run.slim_cmd())
slim -d N={wildcards.N} -d mu={wildcards.mu} -d sh={wildcards.sh} -d rbp={wildcards.rbp} -s {wildcards.seed} -d rep={wildcards.rep} -d "treeseq_file='{output.treeseq_file}'" -d "log_file='{output.log_file}'" -d name=\\"region\\" region.slim
```

These wildcards for each parameter match the free variables in the YAML
configuration file and the SliM command.

### Using SlimFlow in a Snakefile

Then, using these components (the generated targets, target template, and the
SLiM command with wildcards), we can run easily write a Snakefile that runs our
simulations:


```python
from slimflow import GridRuns

if not len(config):
  raise ValueError("config file not specified, use --configfile <config>.yml")

run = GridRuns(config)

# generate all the targets for this YAML file.
all_sims = run.generate_targets()

# a simple, small amount of boilerplate code will run all simulations.
rule slim:
  input: run.script, **run.input
  output: **run.output_template()
  shell: run.slim_cmd()

rule all:
  input: all_sims['filepath']

```

### Working Downstream of Simulation Results

Often, we need to process simulation results in downstream analyses.
`GridRuns.output_template()` can be used to generate the input and output of
Snakefile rules for these downstream analyses. 


