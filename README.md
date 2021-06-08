# envoMatch
Search for peptide isotopic envelopes in ms1 spectra.

# Installation
```bash
git clone https://github.com/weerapana-lab/envoMatch
cd envoMatch
python3 build setup.py
pip install
```

# Usage
```
usage: envoMatch [-h] [--env_co ENV_CO] [--mz_step_margin MZ_STEP_MARGIN]
                 [-t {ms1,mzXML,mzML}] [--ms1_prefix MS1_PREFIX]
                 [-f {input,calculate}] [-s {input,ms1}] [-a ATOM_TABLE]
                 [--plotEnv] [--splitPlots] [--parallel {0,1}]
                 [--nThread NTHREAD] [--overwrite {0,1}] [-v]
                 input_file

positional arguments:
  input_file            ionFinder output file to read.

optional arguments:
  -h, --help            show this help message and exit
  --env_co ENV_CO       Envelope correlation score cutoff.
  --mz_step_margin MZ_STEP_MARGIN
                        Margin above and below envelope in plot.
  -t {ms1,mzXML,mzML}, --file_type {ms1,mzXML,mzML}
                        MS-1 input file type. Default is mzXML.
  --ms1_prefix MS1_PREFIX
                        Append directory to search path for ms1 files. By
                        default only the current working directory is used.
                        Additional directories are searched in the order they
                        are provided.
  -f {input,calculate}, --formula_source {input,calculate}
                        Where should peptide formulas come from? Default is
                        calculate.
  -s {input,ms1}, --pre_scan_src {input,ms1} 
                        Where should precursor scans come from. Chose either
                        the "precursor_scan" column (input) or build precursor
                        list from input MS-1 files (ms1). Default is ms1.
  -a ATOM_TABLE, --atom_table ATOM_TABLE
                        Path to atom table to use in calculating envelopes.
                        Default is: /data/mauraisa/code/envoMatch/db/atom_tabl
                        es/cit_diff_mod_atoms.txt
  --plotEnv             Should plot of envelopes be saved?
  --splitPlots          Split "good" and "bad" envelope plots in separate
                        directories.
  --parallel {0,1}      Chose whether envelope matching should be performed in
                        parallel. Parallel processing is performed on up to
                        the number of logical cores on your system. 1 is the
                        default.
  --nThread NTHREAD     Chose how many threads to use for parllel processing.
                        This option overrides the --parallel option.
  --overwrite {0,1}     Should ionFinder_output be overwritten?
  -v, --verbose         Print verbose output?
  ```
  
  ## `input_file` structure
  
  The input file should be a `.tsv` formated file. At a mininum the input file should have the columns: `parent_file`, `scan`, `sequence`, `formula`, and ` `charge`.
 
 ### Minimal example
