# CoPositivity detection

## Prerequisites:

### - C++17

### - Python3
With numpy and scipy for matrix null space and LP.

### - Install Intel MKL library as a standalone package
1. Go to https://software.seek.intel.com/performance-libraries
2. Register
3. Download Intel MKL
4. tar zxvf l_mkl_2019.0.117.tgz
5. cd folder
6. ./install_GUI.sh
7. Customize installation to install just Intel Math Kernel Library for C/C++ with GNU C/C++ compiler support
8. cd ~/intel/mkl/bin$ && sh ./mklvars.sh intel64

## How to compile:
- **make clean**: Clean build and exec.
- **make**: Compile
- **make total**: Clean, compile and run
- ***IMPORTANT***: 
    - Check your *Intel MKL* and *Python 3* path in the makefile.
    - `.../intel/mkl/lib/intel64` must be in LD_LIBRARY_PATH (make total performs that).
    
## How to run:
Run `run -?` to see help.

```
usage:
  run  options

where options are:
  -n, --dimension <dimension>         Problem dimension (integer) [Required]
  -M, --matrix <matrix_path>          Matrix file (string) [n!=3,5]
  -t <dimacs_t>                       Dimacs t value (integer) [n!=3,5]
  -D, --division <division_method>    Division method (string) ['facet',
                                      'bundfuss', 'zbund'] [Required]
  -?, -h, --help                      display usage information
```

*Remainder:* Matrix dim 3 and dim 5 (Horn) is harcoded, Run them using:
- `run -n 3 -D <division_method>`
- `run -n 5 -D <division_method>`

Division methods:
- ***facet***: Facets and monotonicity (SARTECO18)
- ***bundfuss***: Bundfuss
- ***zbund***: Bundfuss and monotonicity
