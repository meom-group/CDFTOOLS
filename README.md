# CDFTOOLS
  CDFTOOLS is a diagnostic package written in fortran 90 for the analysis of NEMO model output, initialized in  the frame of the DRAKKAR project. It is now available on GitHub under the CeCILL license.

## Using CDFTOOLS

#### Cloning the git repository
To retrieve a copy of the CDFTOOLS source code and create a working directory, run the following on the command line: 

```> git clone https://github.com/meom-group/CDFTOOLS```

#### Compiling CDFTOOLS
There is a Makefile for compiling the sources. The compiler/machines related definitions are supposed to be collected in a `make.macro` file. Some examples of `make.macro` are given in the Macrolib directory and can be used as template for a new compiler or new machine. Then the good practice is to make a link 

```>  ln -sf Macrolib/macro.MACHINE  make.macro```

In the `make.macro` file, the PATH for the netcdf library is specified, as well as compiler name and used flags.  In order to activate netcdf4/HDF5 chunking and deflation ( available in some cdftools), you need to set: 

```NC4=-Dkey_netcdf4 ```

in the make.macro file, otherwise set

```NC4= ```

Then using `make` (or even `make -j n` if you can compile on n cores), you will have the cdftools programs executable available in the bin/ sub directory. The executable files are ignore by git.


#### Running CDFTOOLS
CDFTOOLS is a collection of programs. Every single program performs one or many computation(s) using a set of input files and output the results as a netcdf file, and eventually also gives some results on the standard output. 

CDFTOOLS coding rules imply that there is a build-in documentation foreach cdftool, which is available by just running the tool without any arguments ( or with -h )

## Coding CDFTOOLS
#### Coding rules
##### Syntax
The coding rules are the NEMO coding rules, strictly followed. The idea is that people familiar with NEMO are familiar with CDFTOOLS.
##### Run time behaviour
Any `cdftool`, run without argument or with option -h, should display a short documentation, similar to a unix man page, describing the purpose of the tool, the syntax (arguments,  options, etc...) and giving details on the output files. For some tools, mesh or/and mask files are required to be present in the working directory, with respective name `mesh_hgr.nc`, `mesh_zgr.nc` or `mask.nc` (links are OK). The usage message should indicate the required files.

Example:


```>    cdfcurl```

       usage : cdfcurl -u U-file U-var -v V-file V-var -l levlist [-T] [-8]...
            ... [-surf] [-overf] [-nc4] [-o OUT-file ]
       
      PURPOSE :
        Compute the curl of a vector field, at a specified level.
        If level is specified as 0, assume that the input files are
        forcing files, presumably on A-grid. In this latter case, the
        vector field is interpolated on the C-grid. In any case, the
        curl is computed on the F-point (unless -T option is used).
       
      ARGUMENTS :
        -u U-file U-var : file and variable name for zonal component
        -v V-file V-var : file and variable name for meridional component
        -l levlist    : levels to be processed. If set to 0, assume forcing file
                 in input. Example of recognized syntax :
                   -l "1,10,30"  or -l "1-20" or even -l "1-3,10-20,30-"
                   -l  1 . Note that -l "3-" set a levlist from 3 to the bottom
                   
      OPTIONS :
        -T : compute curl at T point instead of default F-point
        -8 : save in double precision instead of standard simple precision.
        -surf : work with single level C-grid (not forcing)
        -overf : store the ratio curl/f where f is the coriolis parameter
        -nc4 : use netcdf4 output with chunking and deflation 1
        -o OUT-file : specify output file name instead of curl.nc
       
      REQUIRED FILES :
         mesh_hgr.nc
       
      OUTPUT : 
        netcdf file : curl.nc
          variables : socurl or socurlt (if -T option), units : s^-1
             or socurloverf, no units (if -overf option)



