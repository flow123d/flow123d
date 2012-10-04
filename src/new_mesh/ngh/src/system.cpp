#include <cstdio>
#include <cstring>
#include <string>
//#include <stdio.h>
//#include <string.h>
#include  <iostream>

#include "new_mesh/ngh/include/config.h"
#include "new_mesh/ngh/include/system.h"

void mythrow(char *text, int line, const char *fce) {
    static char s[255];
    sprintf(s, "!!!ERROR!!!\nLine: %d, Function: %s\n%s\n", line, fce, text);
    throw s;
}

void warning(char *text, int line, char *fce) {
    static char s[255];
    sprintf(s, "WARNING\nLine: %d, Function: %s\n%s\n", line, fce, text);
    std::cout << s;
}

char *current_directory(char *path) {
    char tmp[MAXPATH];
    getcwd(tmp, MAXPATH);
    strcpy(path, tmp);
    return path;
}

bool GetBool(char *str) {
    if (
            (strcmpi(str, "true") == 0) ||
            (strcmpi(str, "YES") == 0) ||
            (strcmpi(str, "1") == 0)
            ) {
        return true;
    } else {
        return false;
    }
}
