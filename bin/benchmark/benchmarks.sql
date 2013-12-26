# test hosts

CREATE TABLE host (
    id INTEGER PRIMARY KEY,
    hostname TEXT,      # full address or IP
    freqency INTEGER,   # in Hz
    os TEXT,            # win_32, win_64, linux_32, linux_64, possibly other table
);

# one integration test con file, one unit test
# possibly split them
CREATE TABLE test_unit (
    id INTEGER PRIMARY KEY,
    path TEXT,          # path to the file with the test *.con or *.cpp
    size INTEGER,       # size of the problem (integration tests)
);

# individual blocks in the test_unit
CREATE TABLE test_block (
    id INTEGER PRIMARY KEY,
    test_unit_id INTEGER,
    sub_path TEXT,
    leaf INTEGER        # true for leaf blocks
);

# DAG of commits
CREATE TABLE commit (
  id INTEGER PRIMARY KEY,
  parent INTEGER,
  second_parent INTEGER
);

CREATE TABLE test_case (
  id INTEGER PRIMARY KEY,
  commit_id INTEGER,
  host_id INTEGER,
  test_block_id INTEGER,
  
};

CREATE TABLE test_run (
  id INTEGER PRIMARY KEY,
  test_case_id INTEGER,
  
  # setting
  nproc INTEGER,        # number of processors
  #run INTEGER,          # to make data unique
                        # rather update to minimum
                        
  # metrics
  time INTEGER,         # in ms, possibly total of childs
  memory INTEGER,       # in B, possibly total
  calls INTEGER,        # possibly total
  
);