# snake-cnv

Snakemake-based CNV calling pipeline using MrCaNaVar and mrsFAST

## Installation

Get the code.

```bash
git clone --recursive https://github.com/huddlej/snake-cnv.git
cd snake-cnv
```

Install the [latest Python 3 version of Miniconda for your operating system](https://conda.io/miniconda.html).
After installing Miniconda, create a new conda environment with Snakemake installed.

```bash
conda create -n snake-cnv -c bioconda -c conda-forge python=3 conda snakemake pandas
```

If you want to download raw reads using Aspera, download and install the [Aspera command line interface (CLI)](https://downloads.asperasoft.com/en/downloads/62) or find the path to the `ascp` command on your cluster environment.
Edit the `aspera_key_path` parameter in `config.json` to point to the absolute path of your Aspera key.
This is usually something like `<PATH TO ASPERA INSTALLATION>/etc/asperaweb_id_dsa.openssh`.

### Linux-specific installation

Build dependencies for the analysis pipeline.

```bash
make
```

### Mac OS X installation

Install the OpenMP library with [Homebrew](https://brew.sh/) and then build dependencies as follows.

```bash
brew install libomp
make OPENMP_CXX="g++ -lstdc++ -lz -lm -Xpreprocessor -fopenmp -lomp"
```

## Running on a local machine

Run Snakemake on your local machine with conda environment enabled.
The first time you run this will take longer than subsequent runs, since Snakemake will download and install all dependendies in a custom conda environment.

```bash
conda activate snake-cnv
snakemake -w 30 -j 2 --use-conda
```

## Running on a SLURM cluster

```bash
conda activate snake-cnv
snakemake -w 30 -j 2 --use-conda --cluster-config cluster.json --cluster "sbatch --nodes=1 --ntasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time} --job-name='{cluster.name}' --output='{cluster.stdout}' --error='{cluster.stderr}'"
```
