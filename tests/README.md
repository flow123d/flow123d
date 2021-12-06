# Documentation of script `runtest`
Script for running integration tests for application flow123d. Script will
run all tests in given directories (on given yaml files) and will respect
configuration set in `config.yaml` file located in the same directory as
yaml file. (Documentation of the `config.yaml` is specified [below](#documentation-of-configyaml).)


Script `runtest` will use configuration from following sources:
  1. `Global config` - cannot be easily changed, wired in script, default config for all yaml files
  2. file `config.yaml` section `common_config` - (OPTIONAL) easily changable, default config for all yaml files in the same directory
  3. file `config.yaml` section `test_cases` - (OPTIONAL) easily changable, config setting for single test case
  4. parameters for scripts `runtest` such as `-l` ot `-t`
  
Note that configuration from lower level is overriden by higher level configuration. Meaning calling `runtest . -t 30` will set time limit for all cases to 30 seconds.
  
Usage of `runtest` script:
```
runtest [-h] [-a] [-v [VALGRIND]] [--massif] [-p PARALLEL] [--batch]
        [--include [TAG [TAG ...]]] [--exclude [TAG [TAG ...]]]
        [-n CPU] [-q [QUEUE]] [--host HOST] [-t TIME_LIMIT]
        [-m MEMORY_LIMIT] [--root ROOT] [--json JSON] [--list]
        [--dump DUMP] [--log LOG]
        [--random-output-dir [RANDOM_OUTPUT_DIR]] [--no-clean]
        [--no-compare] [--death-test] [--export] [--status-file]
        FOLDER_1 [FOLDER_2 ... FOLDER_N]
        [-- [ARG [ARG ...]]]
```

to list all options run

```sh
$ ./runtest -h
```


# Documentation of `config.yaml`

## Syntax
Yaml file syntax where two section are available:

  1. `common_config` specifying default limits and configuration for all tests in `config.yaml` directory
  2. `test_cases` specifying limits for selected test cases (yaml files).

Both section have almost identical syntax only difference is that section `test_cases` requires key `files` which is array of files for which are next rules applied to.

## Syntax of `common_config`

```yaml
common_config:
  proc: [1, 4]          # int[]     - list of number of MPI processes
                        #   value 0 means no mpi
                        #   to disable test use empty array []
  time_limit:   5.0     # float     - allowed time in seconds
  memory_limit: 800     # int       - allowed memory for ALL processes
  death_test: false     # bool      - true to require the test to fail in
                        #   order to pass the case, default is false
  tags: [foo, bar]      # str[]     - list of tags for the test cases
  args: ['--no-log']    # str[]     - additional arguments passed to the Flow123d binary
  check_rules:          # object    - section where rules can be specified, for every key 
                        #             appropriate module will be loaded 
                        #             (`src/python/scripts/cocmparisons/modules`)
                        #             description of supported comparison rules folows.
      - ndiff:          # object[]  - Perl script `ndiff`, exact match of the text, numerical match of detected 
                        #             numerical valuestext comparison using Perl script 
          files: ["*"]  # str[]     - list of selectors which will select
                        #   what files will be checked using this rule
                        #   by default all files in output directory
                        #   
                        #   To compare only txt files in output dir root:
                        #       files: ["*.txt"]
                        #   To compare only vtu files located in #
                        #   some dir of output dir:
                        #       files: ["*/*.vtu"]
                        #   Wildcards are searched recursively using `Path.walk`, see (`src/python/scripts/prescriptions/__init__.py:AbstractRun._get_ref_output_files).
          r_tol: 1e-3   # float     - optional relative error tolerance
          a_tol: 1e-6   # float     - optional absolute error tolerance
      - regexp:         # Check presence of given regexp in the output file.
          files:        # same as above
          regex: "[Ee]rror"           # matches the word error or Error in line
          substr: "error"             # looks for the word error in a line

          
      - vtkdiff:        # Compute Linf norm of the corresponding datasets using the VTK library. 
                        # allows comparison of binary VTK files.
                        # All datasets from the reference file must be presented. Additional datasets
                        # in the result file are reported as Warning.
          files:        # same as above
          r_tol: 1e-3   # float     - optional relative error tolerance
          a_tol: 1e-6   # float     - optional absolute error tolerance
          interpolate:  # bool      - perform interpolation to the NodeData of the reference mesh before caclulation the norm
              true 

```

## Syntax of `test_cases`

```yaml
# Main list contain one item for every active test.
test_cases:
  - files: flow_implicit.yaml   # str[] | str   - files to which next
                                # configurations is applied to
    
    ## <common_config syntax> ##
```


## Example of `config.yaml`

```yaml
common_config:
  proc: [1]
  time_limit: 30
  memory_limit: 400
  check_rules: 
      - ndiff:
          files: ["*"]
  
test_cases:
  - files:
      - flow_implicit.yaml
      - transport.yaml
    memory_limit: 800
    
  - files: flow_implicit.yaml
    tags: [long-run, problematic]
    proc: [3, 2]
    time_limit: 5.0
    check_rules:
      - ndiff:
            files: ["*/*.vtu"]
            r_tol: 1e-3
            a_tol: 1e-6
      - ndiff:
            files: ["*.vtu"]
            r_tol: 100
            a_tol: 200
    
  - files: flow_implicit.yaml
    tags: [problematic]
    proc: [1]
    time_limit: 5.0
    memory_limit: 800
    check_rules:
      - ndiff:
            files: ["*/*.vtu"]
            r_tol: 1e-3
            a_tol: 1e-6
      - ndiff:
            files: ["*.txt"]
            r_tol: 10
            a_tol: 20
```
