# CDFTOOLS
  CDFTOOLS is a diagnostic package written in fortran 90 for the analysis of NEMO model output, initialized in  the frame of the DRAKKAR project. It is now available on GitHub under the CeCILL license.

## Using CDFTOOLS

#### Cloning the git repository
To retrieve a copy of the SynthPro source code and create a working directory, run the following on the command line: 

```> git clone https://github.com/meom-group/CDFTOOLS```


#### Running CDFTOOLS
CDFTOOLS is a collection of programs. Every single program perform one or many computations using a set of input files and output the results as a netcdf file, and eventually also gives some results on the standard output. 

CDFTOOLS coding rules imply that there is a build-in documentation foreach cdftool, which is available by just running the tool without any arguments ( or with -h )

## Coding CDFTOOLS
#### Coding rules
##### Syntax
The coding rules are the NEMO coding rules, strictly followed. The idea is that people familiar with NEMO are familiar with CDFTOOLS.
##### Run time behaviour
Any cdftool, run without argument or with option -h, should display a short documentation, similar to a unix man page, describing the purpose of the tool, the syntax ( arguments,  options, etc...) and giving details on the output files.

