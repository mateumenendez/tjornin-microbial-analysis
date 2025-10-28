#!/bin/bash

#SBATCH --job-name=Sourcetracking
#SBATCH --partition=compregular
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=10G
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err


sourcetracker2 gibbs -i results/st-biome-class.biom -m results/st-biome-class.map -o results/results-st --sink_rarefaction_depth 0 --source_rarefaction_depth 0 --jobs 32 --per_sink_feature_assignments --restarts 100 --draws_per_restart 5 --diagnostics
