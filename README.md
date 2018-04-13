# omegaMap
Estimating diversifying selection and functional constraint in the presence of recombination

omegaMap is a program for detecting natural selection and recombination in DNA or RNA sequences. It is based on a model of population genetics and molecular evolution. The signature of natural selection is detected using the dN/dS ratio (which measures the relative excess of non-synonymous to synonymous polymorphism) and the signature of recombination is detected from the patterns of linkage disequilibrium. 

The model and the method of estimation are described in

Estimating diversifying selection and functional constraint in the presence of recombination
Daniel J. Wilson and Gilean McVean
Genetics doi:10.1534/genetics.105.044917

This is the source code for Docker Hub image https://hub.docker.com/r/dannywilson/omegamap/

To download and run using Docker:
docker pull dannywilson/omegamap
docker run -u $(id -u):$(id -g) -v $LOCALDIR:/home/ubuntu dannywilson/omegamap $INIFILE

To download and run using Singularity:
singularity pull docker://dannywilson/omegamap
singularity run -H $LOCALDIR:/home/ubuntu omegamap.img $INIFILE

where $LOCALDIR is the local working directory containing the gemma input files and $INIFILE is the configuration file.
