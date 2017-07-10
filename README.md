# snake-cnv

Snakemake-based CNV calling pipeline using MrCaNaVar and mrsFAST

## Installation

Get the code.

```bash
git clone https://github.com/huddlej/snake-cnv.git
```

Build the dependencies.

```bash
make
```

## Running on a SLURM cluster

```bash
snakemake -w 30 -j 2 --cluster-config cluster.json --cluster "sbatch --nodes=1 --ntasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time} --job-name='{cluster.name}' --output='{cluster.stdout}' --error='{cluster.stderr}'"
```
