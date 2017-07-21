# CDFTOOLS
  CDFTOOLS is a diagnostic package written in fortran 90 for the analysis of NEMO model output, initialized in  the frame of the DRAKKAR project (<https://www.drakkar-ocean.eu/>). It is now available on GitHub under the CeCILL license (<http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>).

  CDFTOOLS is an open source package and contributions from other group are welcomed. The Git workflow policy is still to be defined.

  Actual version of CDFTOOLS is known as version 4. (See changes in paragraph *New features in version 4*, below).

## Using CDFTOOLS

#### Cloning the git repository
To retrieve a copy of the CDFTOOLS source code and create a working directory, run the following on the command line: 

```> git clone https://github.com/meom-group/CDFTOOLS ```

#### Compiling CDFTOOLS
All the fortran source are in src/ subdirectory. In src/ there is a Makefile for compiling the sources. The compiler/machines related definitions are supposed to be collected in a `make.macro` file. Some examples of `make.macro` are given in the Macrolib directory and can be used as template for a new compiler or new machine. Then the good practice is to make a link 

```> cd src/```

```> ln -sf ../Macrolib/macro.MACHINE  make.macro```

In the `make.macro` file, the PATH for the netcdf library is specified, as well as compiler name and used flags.  In order to activate netcdf4/HDF5 chunking and deflation ( available in some cdftools), you need to set: 

```NC4=-Dkey_netcdf4 ```

in the make.macro file, otherwise set

```NC4= ```

In order to activate CMIP6 variable naming convention (for input files), you need to set:

```CMIP6=-Dkey_CMIP6 ```

Then using `make` (or even `make -j n` if you can compile on n cores), you will have the cdftools programs executable available in the bin/ sub directory. The executable files are ignore by git.


#### Running CDFTOOLS
CDFTOOLS is a collection of programs. Every single program performs one or many computation(s) using a set of input files and output the results as a netcdf file, and eventually also gives some results on the standard output. 

CDFTOOLS coding rules imply that a `usage message` is displayed when just running the tool without any arguments ( or with -h ). At the moment it is the only up to date documentation. 

As CDFTOOLS is a collection of programs, a full diagnostic of model output can be performed writing a script using a sequence of tools. This is done for example in the Drakkar Monitoring Tools (DMONTOOLS, soon available on GitHub!).

## Coding CDFTOOLS
#### Coding rules
##### Syntax
The coding rules are the NEMO coding rules, strictly followed. The idea is that people familiar with NEMO are familiar with CDFTOOLS. In DEV_TOOLS/ some template fortran line are available for program, modules, routine or function headers. Also a template for the `usage message`.
##### Run time behaviour
Any `cdftool`, run without argument or with option -h, should display a short documentation (`usage message`), similar to a unix man page, describing the purpose of the tool, the syntax (arguments,  options, etc...) and giving details on the output files. For some tools, mesh or/and mask files are required to be present in the working directory, with respective name `mesh_hgr.nc`, `mesh_zgr.nc` or `mask.nc` (links are OK). The usage message should indicate the required files.

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

##### Improving/modifying existing tool
 It is possible to improve (of course !) or modify any tools, but <u>one important law to respect</u> is that the modified tool should still be able to be used with previous syntax, in order to avoid breaking of existing scripts using CDFTOOLS. If for some reason, this is not possible, then a discussion must be done to reach to a common decision. Eventually, some old options must be documented as osbolete in the usage message, which means that they may be removed from a future release. 

## New features in version 4
#### Modified user interface
 * All arguments are passed with a `-key` switch. No more `free` arguments. Example `cdfmoy -l fich1.nc fich2.nc`
 * Add `-o` and `-nc4` options in all tools (when relevant). With `-o` the default output name can be changed, allowing easier paralellisation. With `-nc4` output file will use NetCdf4/Hdf5 format with chunking and deflation level 1.
 * Use of environment variable CDFT_xxx for overriding  the default names of auxiliary files such as mesh_hgr.nc, mask.nc etc...so far there is support for 

   CDFT_MESH_HGR
   CDFT_MESH_ZGR
   CDFT_MASK
   CDFT_BASINS
   CDFT_COORD

#### Support for vvl simulations
 * When relevant, the switch `-vvl` indicates that the vertical metrics is time-varying. Therefore, CDFTOOLS assume that the vertical metrics is saved in the same file than the data.

#### Support for CMIP6 naming convention
 * When the code is compiled with CPP key key_CMIP6 set, the default variable names are taken form modcdfnames_CMIP6.h90 instead of the standard DRAKKAR names.

#### Simplification
 * The codes have been cleaned for obsolescences. Coding rules were reinforced.
 * Obsolete tools were removed or merged as options in more generic tools. 

#### Improved documentation.
 * Gathering the `usage` message into man pages still works (`make man`). Readibility of the man pages is improved by grouping the tools by category. The `usage` messages have been reviewed in order to give better information.
 * The man pages are automaticaly translated to an html document that can be vizualized from any browser.

#### Back to release 3:
 * The last version 3's release has been tagged as v3.0.2. Use this tag if you want to stay at version 3.

