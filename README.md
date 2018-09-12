# snake-cnv

Snakemake-based CNV calling pipeline using MrCaNaVar and mrsFAST

## Installation

Get the code.

```bash
git clone --recursive https://github.com/huddlej/snake-cnv.git
```

Build dependencies for the analysis pipeline.

```bash
make
```

If you already have Python 3 installed, you only need to install Snakemake.

```bash
pip3 install snakemake
```

If you do not have Python 3, install the [latest version of Miniconda for Python 3 for your machine](https://conda.io/miniconda.html).
After installing Miniconda, install Snakemake with conda.

```bash
conda install -c bioconda snakemake
```

## Running on a local machine

Run Snakemake on your local machine with conda environment enabled.
The first time you run this will take longer than subsequent runs, since Snakemake will download and install all dependendies in a custom conda environment.

```bash
snakemake -w 30 -j 2 --use-conda
```

## Running on a SLURM cluster

```bash
snakemake -w 30 -j 2 --use-conda --cluster-config cluster.json --cluster "sbatch --nodes=1 --ntasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time} --job-name='{cluster.name}' --output='{cluster.stdout}' --error='{cluster.stderr}'"
```
