#ifndef problemH
#define problemH

#include "config.h"
#include <stdio.h>
#include "system.h"

typedef enum ComputeTypes {
    ct_none,
    ct_ratio1
} TComputeType;

#define _VERSION_ "0.1.1"

extern double epsilon;

class TProblem {
private:
    char ini_fname[MAXPATH];
    char mesh_fname[MAXPATH];
    char output_fname[MAXPATH];

    bool pause_after_run;

    bool compute_bb;
    bool compute_vb;
    bool compute_vv;

    TComputeType ct_11;
    TComputeType ct_12;
    TComputeType ct_13;

    TComputeType ct_22;
    TComputeType ct_23;

    void Initialize();
    void ReadFile();

public:
    TProblem(char*);
    ~TProblem();

    bool getBB();
    bool getVB();
    bool getVV();

    TComputeType getCt_11();
    TComputeType getCt_12();
    TComputeType getCt_13();

    TComputeType getCt_22();
    TComputeType getCt_23();

    bool getPause_after_run();

    char* getMesh_fname();
    char* getOutput_fname();
};

TComputeType GetCT(const char*);

#endif
