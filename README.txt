DISCLAIMER
-------------------------------------------------------------------------------
MLML has become part of METHPIPE (https://github.com/smithlabcode/methpipe), 
and is actively managed in METHPIPE. This stand-alone MLML repository is no 
longer maintained as of Feb 2018. 

For more details, please read the user manual in docs/manual.pdf

PREREQUISITES
-------------------------------------------------------------------------------
GNU Compiler Collection (GCC)
  This software requires a C++ complier to build. The most recent version of
  GCC is highly recommended. You can find it at http://gcc.gnu.org/

GNU Scientific Library (GSL)
  This software package requires a functioning installation of the GNU 
  Scientific Library (GSL). If you don't already have this installed, you 
  will need to download and install it from http://www.gnu.org/software/gsl/

INSTALLATION
-------------------------------------------------------------------------------
This software package is developed in UNIX environment. It has been tested
in Mac OSX and Linux.

  (1) To build binaries, type the following command in mlml root directory:
  make

  (2) To install binaries, after building type the following command:
  make install

If no ERROR messages prompts during the two steps, then your copy of mlml
has been installed successfully.

USAGE
-------------------------------------------------------------------------------
-o, -output
  Specify the output file path. If this option is not given, then the output
  will be directed to screen output.

-u, -bsseq
  Specify the input file path for BS-seq file. The input file format should be
  BED format.

-h, -tabseq
  Specify the input file path for Tab-seq file. The input file format should be
  BED format.

-m, -oxbsseq
  Specify the input file path for oxBS-seq file. The input file format should be
  BED format.

* Note: at least two of -u, -h and -m options are required.

-t, -tolerance
  Specify the convergence tolerance of EM algorithm. Default value is 1e-10.
  Larger value will yield more accurate result, with the cost of longer runtime.

-a, -alpha
  Specify the alpha of binomial test confidence interval. Default value is 0.05.
  Smaller value will make the test more stringent.

-v, -verbose
  Turn on verbose mode. Default off. If turned on, the program will print
  detailed information to the screen of run statistics. Recommended.

EXAMPLES
-------------------------------------------------------------------------------
(1) all 3 types of inputs are available:

  mlml -v -u FILE_BS -m FILE_OX -h FILE_TAB -o OUTPUT_FILE

  This will make the program load methylation data from FILE_BS, FILE_OX and
  FILE_TAB respectively, and write output to OUTPUT_FILE

(2) only BS-seq data and Tab-seq data are available:

  mlml -v -u FILE_BS -h FILE_TAB -o OUTPUT_FILE

  The program will load methylation data from FILE_BS and FILE_TAB.

INPUT FILE FORMAT
-------------------------------------------------------------------------------
All input files for -u, -m and -h options should be in BED format. Detailed
file format is discussed in manual. Followed is an example:

chr1	3001345 3001346 CpG:9   0.777777777778  +

Here the 4th column indicates that the genomic coordinate is a CpG site, and
the number of reads covering this site is 9. The 5th column is methylation
level of the site.

CONTACTS
-------------------------------------------------------------------------------
Jenny Qu
jqu@usc.edu

Meng Zhou
mengzhou@usc.edu

Andrew D. Smith
andrewds@usc.edu

Copyright and License Information
-------------------------------------------------------------------------------
Copyright (C) 2013
University of Southern California,
Jenny Qu, Meng Zhou and Andrew D. Smith
  
Authors: Jenny Qu, Meng Zhou and Andrew D. Smith
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
