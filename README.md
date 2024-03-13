# SlimFlow — A lightweight framework for running SLiM simulations with Snakemake

SlimFlow is a small Python module that makes running and processing [SLiM
simulations](https://messerlab.org/slim/) with
[Snakemake](https://snakemake.readthedocs.io/en/stable/) much easier. It will
also work with [msprime](https://tskit.dev/software/msprime.html), though this
is currently in development. I've used this module for several projects (e.g.
Buffalo and Coop 2020, Buffalo and Kern 2024), and split it out so others can
use. 

## Installation

SlimFlow is not yet on PyPi but in the meantime, you can use:

```
$ pip install git+https://github.com/vsbuffalo/slimflow
```

## Example

To step through a full example, see the documentation. Here is a brief look 
at the API, so you can get some idea of how it can be used in a Snakefile:

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

Then, the targets files (e.g. simulation results or their downstream processed
files) can be generated using `GridRuns.generate_targets()`. This produces a
[Polars](https://pola.rs) dataframe of all simulation results, created by
creating all parameter combinations for all replicates, seeding them
appropriately.


```python
>>> run.generate_targets().head(2)
shape: (2, 10)
┌───────────────────────────────────┬──────┬───────────┬──────┬───┬───────────────────────────────────┬─────┬─────────────────────┬──────────────┐
│ key                               ┆ N    ┆ mu        ┆ sh   ┆ … ┆ filename                          ┆ rep ┆ seed                ┆ suffix       │
│ ---                               ┆ ---  ┆ ---       ┆ ---  ┆   ┆ ---                               ┆ --- ┆ ---                 ┆ ---          │
│ str                               ┆ i64  ┆ f64       ┆ f64  ┆   ┆ str                               ┆ i64 ┆ i64                 ┆ str          │
╞═══════════════════════════════════╪══════╪═══════════╪══════╪═══╪═══════════════════════════════════╪═════╪═════════════════════╪══════════════╡
│ runs___bgs___N__1000___mu__2e-08… ┆ 1000 ┆ 2.0000e-8 ┆ 0.01 ┆ … ┆ rep_0__seed_7138484576005690179_… ┆ 0   ┆ 7138484576005690179 ┆ treeseq.tree │
│ runs___bgs___N__1000___mu__2e-08… ┆ 1000 ┆ 2.0000e-8 ┆ 0.01 ┆ … ┆ rep_0__seed_7138484576005690179_… ┆ 0   ┆ 7138484576005690179 ┆ log.tsv.gz   │
└───────────────────────────────────┴──────┴───────────┴──────┴───┴───────────────────────────────────┴─────┴─────────────────────┴──────────────┘
```

We can save these targets, and take a look at how SlimFlow organizes simulation results
in a tidy way:

```python
>>> targets = run.generate_targets()

# targets were generated for both suffices for this simulation replicate:
>>> targets['filepath'][0]
'runs/bgs/N__1000/mu__2e-08/sh__0.01/rbp__1e-08/rep_0__seed_7138484576005690179__treeseq.tree'

```

Then, the `filepath` column has the relevant information on each file produced
by the simulation. This Polars dataframe can then be joined in with downstream
results.

SlimFlow removes a ton of boilerplate code, allowing one to run and manage SLiM
results effortlessly. This is a small Snakefile that shows how to generate all
of these simulations:


```python
from slimflow import GridRuns

if not len(config):
    raise ValueError("config file not specified, use --configfile <config>.yml")

run = GridRuns(config)

# Generate all the targets for this YAML file.
sims_df = run.generate_targets()

# The SLiM rule: a simple, small amount of boilerplate code will run all
# simulations.
rule slim:
  input: run.script, **run.input
  output: **run.target_template()
  shell: run.slim_cmd()

rule all:
  input: sims_df['filepath']
```

## Similar Projects

 - [stdpopsim](https://github.com/popsim-consortium/stdpopsim) is a larger framework
   for running simulations for a variety of different species. This is recommended if
   you need specific demographic models for species histories' etc.

 - [slim_magic](https://github.com/andrewkern/slim_magic) is a neat bit of iPython
   magic that makes working with SLiM simulations in Jupyter notebooks easy.
