# OnlineDNA UDP workflow
Demultiplexes and maps Nanopore reads barcoded with Illumina UDP-style adapters

## Software requirements and dependencies
Run the `install_singularity.sh` script to install singularity either system-wide or only for the current user. Everything is tailored (and tested) for a 64-bit x86 system running Ubuntu Linux 18.04 or later, or Debian Buster or later.

## Building singularity image
In the repository root `/` hit `make base` and afterwards `make main`. This will build the base image with all software and dependencies and on top of that a main image with also scripts included (mainly to avoid rebuilding image repeatedly during development).

## Running container and workflow
Run `singularity run singularity/onlinedna.sif`. 

For an interactive console start an instance with `singularity instance start singularity/onlinedna.sif onlinedna` and `singularity shell instance://onlinedna`. 

Note: By default singularity makes the `/tmp` folder, home folder, and the current WD available from the host inside the container with the same paths.

## Development tips
To reduce waiting time for rebuilding image:
 - Unlike Docker, singularity does not cache individual steps, so test installation of new software by adding to the %post section in the `main.singularity` image file and rebuild with `make`. When succesful move it to the `base.singularity` image file, then rename the `.sif` file to `base.sif` or rebuild base again.
 - For installing additional R packages use `renv::install(pkg)`, then `renv::snapshot()` and then rebuild base image to include them.
 - When developing scripts just live mount the scripts folder (`--bind scripts/:/opt/scripts/`) when launching the container instead of rebuilding image repeatedly.
