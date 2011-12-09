#include "config.h"
#include <iostream>
#include "simpleini.h"
#include "problem.h"
#include "system.h"

void TProblem::Initialize() {
    strcpy(ini_fname, "");
    strcpy(mesh_fname, "");
    strcpy(output_fname, "");

    compute_bb = false;
    compute_vb = false;
    compute_vv = false;

    ct_11 = ct_none;
    ct_12 = ct_none;
    ct_13 = ct_none;
    ct_22 = ct_none;
    ct_23 = ct_none;

    pause_after_run = true;
}

TProblem::TProblem(char* ini) {
    Initialize();

    std::cout << "Opening input file... ";

    FILE *F = fopen(ini, "rt");
    if (F == NULL) {
        mythrow((char*) "Couldn't open input ini file.", __LINE__, __FUNC__);
    }
    strcpy(ini_fname, ini);
    fclose(F);

    std::cout << "OK.\n";
    ReadFile();
}

void TProblem::ReadFile() {
    std::cout << "Reading input file... ";

    CSimpleIni ini; // no UTF8, no Multi Key/Lines
    if (ini.LoadFile(ini_fname) < 0) {
        mythrow((char*) "Failed to open the ini file.", __LINE__, __FUNC__);
    }

    // RUN
    pause_after_run = ini.GetBoolValue("Run", "pause_after_run", "true");

    // FILES
    strcpy(mesh_fname, ini.GetValue("Files", "mesh", ""));
    strcpy(output_fname, ini.GetValue("Files", "output", ""));

    // COMPUTATION
    compute_bb = ini.GetBoolValue("Computation", "bb", "false");
    compute_vb = ini.GetBoolValue("Computation", "vb", "false");
    compute_vv = ini.GetBoolValue("Computation", "vv", "false");

    ct_11 = GetCT(ini.GetValue("Computation", "ct_11", "none"));
    ct_12 = GetCT(ini.GetValue("Computation", "ct_12", "none"));
    ct_13 = GetCT(ini.GetValue("Computation", "ct_13", "none"));
    ct_22 = GetCT(ini.GetValue("Computation", "ct_22", "none"));
    ct_23 = GetCT(ini.GetValue("Computation", "ct_23", "none"));

    // CONSTANTS
    epsilon = atof(ini.GetValue("Constants", "epsilon", "1e-8"));

    std::cout << "OK\n";
}

TComputeType GetCT(const char* str) {
    if (strcmpi(str, "none") == 0) return ct_none;
    if (strcmpi(str, "ratio1") == 0) return ct_ratio1;

    return ct_none;
}

TProblem::~TProblem() {
    //  delete mesh;
}

bool TProblem::getBB() {
    return compute_bb;
}

bool TProblem::getVB() {
    return compute_vb;
}

bool TProblem::getVV() {
    return compute_vv;
}

TComputeType TProblem::getCt_11() {
    return ct_11;
}

TComputeType TProblem::getCt_12() {
    return ct_12;
}

TComputeType TProblem::getCt_13() {
    return ct_13;
}

TComputeType TProblem::getCt_22() {
    return ct_22;
}

TComputeType TProblem::getCt_23() {
    return ct_23;
}

bool TProblem::getPause_after_run() {
    return pause_after_run;
}

char* TProblem::getMesh_fname() {
    return mesh_fname;
}

char* TProblem::getOutput_fname() {
    return output_fname;
}
