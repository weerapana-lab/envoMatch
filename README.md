# envoMatch
Search for peptide isotopic envelopes in MS1 spectra.

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
                        Default is: envoMatch/db/atom_tables/cit_diff_mod_atoms.txt
  --plotEnv             Should plot of envelopes be saved?
  --splitPlots          Split "good" and "bad" envelope plots in separate
                        directories.
  --parallel {0,1}      Chose whether envelope matching should be performed in
                        parallel. Parallel processing is performed on up to
                        the number of logical cores on your system. 1 is the
                        default.
  --nThread NTHREAD     Chose how many threads to use for parallel processing.
                        This option overrides the --parallel option.
  --overwrite {0,1}     Should ionFinder_output be overwritten?
  -v, --verbose         Print verbose output?
```
  
## input_file structure
  
The input file should be a `.tsv` formatted file. At a minimum the input file should have the columns: `parent_file`, `scan`, `sequence`, `formula`, and `charge`.
 
### Minimal example input file

| parent_file | scan | sequence | formula | charge |
| ----------- | ---- | -------- | ------- | ------ |
| file.mzML | 47677 | AEMMELNDR*FASYIEK | C84H130N20O29S2 | 2 |
| file.mzML | 7415  | AGDKAGR*AGAGmPPYHR  | C72H113N25O23S |  4 |
| file.mzML | 41163 | AGFLVTAR*GGSGIVVAR  | C72H122N22O21 | 2 |
| file.mzML | 8204  | AGR*AGAGmPPYHR  | C57H88N20O17S | 3 |
| file.mzML | 28728 | ALAR*EVDLKDYEDQQK | C82H132N22O31 | 2 |
| file.mzML | 47838 | ALQEqLATqR*EAIILAR  | C83H143N23O29 | 2 |

### output
`envoMatch` Adds 2 columns to the input file and writes the file `<input_name>_env.tsv`

1. `good_envelope`: A boolean representing whether the precursor envelope passed all the quality tests.
2. `env_score`: The similarity score comparing the theoretical envelope to the observed.

## Example

To run the `envoMatch` example:
```bash
cd examples
envoMatch input.tsv
```
 The expected output is `examples/output/input_env.tsv`

