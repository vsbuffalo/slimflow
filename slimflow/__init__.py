import sys
import itertools
import os
from os.path import join
import numpy as np
import polars as pl


# this is to guarantee a random seed between 0 and sys.maxize
# does not overflow.
MAX_REPLICATES = 100_000


class GridRuns:
    """
    `GridRuns` is a class to aid in setting up multiple simulations
    across a grid (Cartesian product) of parameter combinations. It will
    translate a YAML file like this:

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

    Into an object containing this information, that can be used to generate
    SLiM targets, etc. Let's look at the loading a YAML config file into a
    set of simulation parameters:

    ```python
    >>> import yaml
    >>> from slimflow import GridRuns
    >>> config = yaml.safe_load(open('tests/data/test_config.yml'))
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

    >>> targets = run.generate_targets()

    >>> targets.columns
    ['key', 'N', 'mu', 'sh', 'rbp', 'filepath', 'filename', 'rep', \
'seed', 'suffix']

    # targets were generated for both suffices for this simulation replicate:
    >>> targets['filepath'][0]
    'runs/bgs/N__1000/mu__2e-08/sh__0.01/rbp__1e-08/rep_0__seed_7138484576005612783__treeseq.tree'
    >>> targets['filepath'][1]
    'runs/bgs/N__1000/mu__2e-08/sh__0.01/rbp__1e-08/rep_0__seed_7138484576005612783__log.tsv.gz'

    # slimflow also generates the proper SLiM command line call to pass
    # these parameters in:
    >>> print(run.slim_cmd())  # doctest: +ELLIPSIS
    slim -d N={wildcards.N} -d mu={wildcards.mu} -d sh={wildcards.sh}...bgs.slim

    >>> print(run.target_template())  # doctest: +ELLIPSIS
    {'treeseq_file': 'runs/bgs/N__{N}/mu__{mu}/sh__{sh}/...}

    >>> print(run.target_template(outdir='processed', suffices={'tsv': '.tsv'},\
                                  filename='results'))
    {'tsv': 'runs/bgs/processed/N__{N}/mu__{mu}/sh__{sh}/.../results.tsv'}

    ```


    This class uses a fairly general YAML format to store simulations in a way
    where they can be quickly found if needed. The format used is:

    This is the parameter-to-directory scheme, and it puts some mild
    constraints on what can be a free parameter (e.g., we can't have file paths
    become the name of a subdirectory).

    Args:
        config (dict): The configuration dictionary for the SLiM runs.

    Attributes:
        config (dict): The configuration dictionary for the SLiM runs.

        name (str): The main run name, which should be unique to the fixed \
                parameters.

        script (str): The path to the SLiM script.

        seed (int): The random seed value.

        rng (numpy.random.Generator): The random number generator object.

        dir (str): The base directory to store the simulations.

        suffices (dict): A dictionary of suffices for the output files.

        fixed (dict): The parameters that do not change across the simulations.

        input (dict): The input files that are also fixed but will trigger a \
          re-run.

        variable (dict): The parameters that vary across simulations.

        nreps (int): The number of repetitions for each parameter combination.

        basedir (str): The base directory path for the simulations.

        simdir (str): The directory path for storing the simulation results.

    """

    def __repr__(self):
        if self.nreps >= MAX_REPLICATES:
            msg = f"nreps = {self.nreps}, but nreps must be < {MAX_REPLICATES}"
            raise ValueError(msg)

        # Total number of combinations
        total_combinations = np.prod([len(v) for v in self.variables.values()])

        # Variable summaries, each on its own line
        var_summaries = [
            f"  {k} ∈ {{{', '.join(map(str, v))}}}"
            for k, v in self.variables.items()
        ]
        var_summary_str = "\n".join(var_summaries)
        total_sims = total_combinations * self.nreps
        return (
            f"GridRuns({self.name}):\n"
            f"Number of parameter combinations (𝑘): {total_combinations}\n"
            f"Number of replicates (𝑛): {self.nreps}\n"
            f"Total simulations (𝑘 × 𝑛): {total_sims}\n"
            f"Variables:\n{var_summary_str}"
        )

    def __init__(self, config):
        self.config = config
        self.name = config["name"]
        self.script = config["script"]
        self.seed = int(config["seed"])
        self.rng = np.random.default_rng(self.seed)
        self.dir = config["dir"]
        self.suffices = config["suffices"]
        self.fixed = config.get("fixed", {})
        self.input = config.get("input", {})
        self.variables = config["variables"]
        msg = "'rep' cannot be used as a variables parameter"
        assert "rep" not in self.variables, msg
        self.nreps = config["nreps"]
        self.basedir = join(self.dir, self.name)

        # this is the dataframe of targets, once generated
        self.targets = None

        simdir = config.get("simdir", None)
        if simdir is not None:
            self.simdir = join(self.basedir, simdir)
        else:
            self.simdir = self.basedir

    def generate_targets(
        self,
        suffices=None,
        filename=None,
        outdir=None,
        nreps=None,
        create=True,
    ):
        """
        Generate all the target filenames for the SLiM runs.

        Args:
            suffices (dict, optional): A dictionary of suffices for the output\
             files. If not provided, the suffices from the config will be used.

            filename (str, optional): A filename for a particular parameter \
            combinations. If `None` (the defualt), the filename will be the \
            automatically-generated filename containing the seed and replicate\
            number. Specifying the filename will drop these per-replicate bits\
            of information, which is useful when multiple replicates are to be\
            processed together and the results output to a single file.

            outdir (str, optional): The output directory path. If not \
                    provided, the simulation directory will be used.

            nreps (int, optional): The number of repetitions for each \
                    parameter combination. If not provided, the value \
                    from the config will be used.

            create (bool, optional): Whether to create the directories if \
                    they don't exist. Defaults to True.

        Returns:
            A list of target filenames.
        """
        # initialize the rng
        rng = np.random.default_rng(self.seed)

        # make all grid combinations of variable parameters
        all_run_params = param_grid(self.variables)
        if suffices is not None:
            assert isinstance(suffices, dict), "suffices must be a dict"
        else:
            suffices = self.suffices
        if outdir is None:
            # use the simulation directory
            outdir = self.simdir
        else:
            # use this directory
            outdir = join(self.basedir, outdir)

        # allow custom nreps
        nreps = self.nreps if nreps is None else nreps

        targets = []
        targets_df = []
        for params in all_run_params:
            # generate the directory structure on the
            # variable params
            dir = params_to_dirs(params, basedir=outdir, create=create)
            # get a random seed
            seed = random_seed(rng)

            row = dict(key=dir.replace("/", "___"))
            row = row | params
            for rep in range(nreps):
                seed += rep  # this way, adding reps won't mess up seeds
                for name, suffix in suffices.items():
                    if filename is not None:
                        output_filename = f"{filename}{suffix}"
                    else:
                        output_filename = filename_pattern(
                            suffix, rep, seed, None
                        )

                    filepath = join(dir, output_filename)
                    targets.append(filepath)

                    suffix_rep_row = row | dict(
                        filepath=filepath,
                        filename=output_filename,
                        rep=rep,
                        seed=seed,
                        suffix=suffix,
                    )
                    targets_df.append(suffix_rep_row)

        targets_df = pl.DataFrame(targets_df)
        return targets_df

    def slim_cmd(self):
        """
        Generate the SLiM command for running the simulations.

        Fixed parameters, input files, and the name are passed directly to the
        SLiM call here and are not part of the wildcards.

        Returns:
            The SLiM command.
        """
        # all the same fixed parameters are passed manually
        manual = dict(**self.fixed, **self.input, name=self.name)

        # we get one sample run, to extract the types
        all_run_params = list(param_grid(self.variables))[0]

        output = self.suffices
        return slim_call(
            all_run_params,
            self.script,
            add_seed=True,
            add_rep=True,
            output=output,
            manual=manual,
        )

    def grouped_targets(self, outdir=None, suffices=None):
        """
        Create a wildcard function for Snakemake that groups targets
        by their parameters, collecting all replicates into a group.
        This allows for easier downstream processing.
        """
        targets_df = self.generate_targets(suffices=suffices, outdir=outdir)

        def input_filenames(wildcards):
            wildcard_dict = dict(wildcards)

            # Remove non-parameter wildcards
            wildcard_dict.pop("rep", None)
            wildcard_dict.pop("seed", None)

            # Convert wildcard values to the appropriate types
            for param, value in wildcard_dict.items():
                if param in self.variables:
                    param_type = type(self.variables[param][0])
                    wildcard_dict[param] = param_type(value)

            # Filter the targets dataframe based on the wildcard values
            filtered_df = pl.DataFrame(targets_df)
            for param, value in wildcard_dict.items():
                filtered_df = filtered_df.filter(pl.col(param) == value)

            # Get the list of input files
            input_files = filtered_df["filepath"].to_list()

            return input_files

        return input_filenames

    def target_template(
        self,
        outdir=None,
        suffices=None,
        fix_params=False,
        filename=None,
        pin_params=None,
    ):
        """

        Create a dictionary of expected output files for all the suffices, with
        wildcards. This is designed to be used directly in the `output` entry
        of a Snakemake rule (which takes can take a set of key value
        arguments), as `output: **run.target_template()`.

        By default, this will produce outputs for the suffixes in the config in
        the base directory set by the config. However, sometimes simulation
        results are needed for further downstream results; in this case, set
        outdir and suffices to something other than None.

        Args:

            outdir (str, optional): The output directory path. If not \
              provided, the simulation directory will be used.

            suffices (dict, optional): A dictionary of suffices for the \
                    output files. If not provided, the suffices from the \
                    config will be used.

            filename (str, optional): A filename for a particular \
                    parameter combinations. If `None` (the defualt), the \
                    filename will be the automatically-generated filename \
                    containing the seed and replicate number. Specifying \
                    the filename will drop these per-replicate bits of \
                    information, which is useful when multiple replicates \
                    are to be processed together and the results output to \
                    a single file.

            pin_params (dict, optional): A dictionary of parameter values to \
              pin. If provided, these parameters will be fixed to the \
              specified values. This is useful if downstream analyses only \
              should process a subset of parameters.

        Returns:
            A dictionary of expected output filenames with wildcards.
        """
        param_wildcards = {}

        if pin_params is not None:
            invalid_pin_params = [
                p for p in pin_params if p not in self.variables
            ]
            if len(invalid_pin_params):
                invalid_str = ", ".join([f"'{x}'" for x in invalid_pin_params])
                raise ValueError(
                    f"parameters {invalid_str} are not in the variables"
                )

        for param_name in self.variables:
            if pin_params is not None and param_name in pin_params:
                param_wildcards[param_name] = pin_params[param_name]
            else:
                if fix_params:
                    param_name = "{" + param_name + ""
                param_wildcards[param_name] = "{" + param_name + "}"

        outputs = {}
        if suffices is not None:
            assert isinstance(suffices, dict), "suffices must be a dict"
        else:
            suffices = self.suffices
        if outdir is None:
            # use the simulation directory
            outdir = self.simdir
        else:
            # use this directory
            outdir = join(self.basedir, outdir)

        for name, suffix in suffices.items():
            if filename is None:
                output_filename = filename_pattern(
                    suffix, "{rep}", "{seed}", None
                )
            else:
                output_filename = f"{filename}{suffix}"
            dir = params_to_dirs(param_wildcards, basedir=outdir)
            outputs[name] = join(dir, output_filename)
        return outputs


def param_grid(params):
    """
    Generate a Cartesian product parameter grid from a dict of grids per
    parameter.

    Args:
        params (dict): A dictionary where keys are parameter names and values
                       are lists of parameter values.

    Returns:
        A list of dictionaries, each representing a combination of parameters.
    """
    grid = []
    for param, values in params.items():
        if len(values):
            grid.append([(param, v) for v in values])
        else:
            grid.append([(param, "")])

    def convert_to_dict(x):
        # and package seed if needed
        x = dict(x)
        return x

    return map(convert_to_dict, itertools.product(*grid))


def random_seed(rng=None):
    """
    Generate a random seed value.

    Args:
        rng (numpy.random.Generator, optional): A random number generator
             object. If not provided, numpy.random will be used.

    Returns:
        int: A random seed value.
    """
    seed_max = sys.maxsize - MAX_REPLICATES
    if rng is None:
        return np.random.randint(0, seed_max)
    return rng.integers(0, seed_max)


def filename_pattern(suffix, rep, seed, params=None):
    """
    Build a filename from parts.

    Args:
        suffix (str): The suffix of the filename.
        rep (int): The repetition number.
        seed (int): The random seed value.
        params (dict, optional): A dictionary of variable parameters. If
         provided, these parameters will be used to construct a more unique ID.

    Returns:
        str: The constructed filename pattern.
    """
    if params is not None:
        param_str = "_".join([v + "{" + v + "}" for v in params]) + "_"
    else:
        param_str = ""
    pattern = f"{param_str}rep_{rep}__seed_{seed}__{suffix}"
    return pattern


def slim_call(
    params,
    script,
    slim_cmd="slim",
    output=None,
    add_seed=True,
    add_rep=True,
    manual=None,
):
    """
    Create a SLiM call prototype for Snakemake.

    Args:
        params (dict): A dictionary of sample parameters, so types can be
          inferred.
        script (str): The path to the SLiM script.
        slim_cmd (str, optional): The path to the SLiM executable. Defaults to
          "slim".
        output (dict, optional): A dictionary of output names and values.
        add_seed (bool, optional): Whether to pass in the seed with '-s
          <seed>'. Defaults to True.
        add_rep (bool, optional): Whether to pass in the repetition number with
          '-d rep={wildcards.rep}'. Defaults to True.
        manual (dict, optional): A dictionary of manual items to pass in.

    Returns:
        str: The SLiM call command.
    """
    param_types = {k: type(v) for k, v in params.items()}
    call_args = []
    for p, val_type in param_types.items():
        is_str = val_type is str
        if is_str:
            val = f'\\"{{wildcards.{p}}}\\"'
        else:
            val = f"{{wildcards.{p}}}"
        call_args.append(f"-d {p}={val}")
    add_on = ""
    if manual is not None:
        add_on = []
        for key, val in manual.items():
            if isinstance(val, str):
                add_on.append(f'-d {key}=\\"{val}\\"')
            else:
                add_on.append(f"-d {key}={val}")
        add_on = " " + " ".join(add_on)
    if add_seed:
        call_args.append("-s {wildcards.seed}")
    if add_rep:
        call_args.append("-d rep={wildcards.rep}")
    if output is not None:
        for name, file in output.items():
            call_args.append(f"-d \"{name}='{{output.{name}}}'\"")
    full_call = f"{slim_cmd} " + " ".join(call_args) + add_on + " " + script
    return full_call


def create_directories(path):
    """
    Create directories if they don't exist.

    Args:
        path (str): The directory path to create.
    """
    if not os.path.exists(path):
        os.makedirs(path)


def params_to_dirs(params, basedir=None, create=False):
    """
    Convert parameters to directory paths.

    Args:
        params (dict): A dictionary of parameters.

        basedir (str, optional): The base directory path. If provided, the
          generated directory will be appended to it.

        create (bool, optional): Whether to create the directories if they
           don't exist. Defaults to False.

    Returns:
        str: The generated directory path.
    """
    dir = join(*[f"{k}__{v}" for k, v in params.items()])
    if basedir is not None:
        dir = join(basedir, dir)
    if create:
        create_directories(dir)
    return dir
