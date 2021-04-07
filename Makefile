all: main
main:
	cd singularity/ && singularity build --fakeroot -F onlinedna.sif main.singularity

base:
	cd singularity/ && singularity build --fakeroot -F base.sif base.singularity

clean:
	cd singularity/ && rm -rf base.sif onlinedna.sif build
