
#include <petsc.h>
#include "system/system.hh"
#include "io_namehandler.hh"

IONameHandler* IONameHandler::instance = NULL;

IONameHandler* IONameHandler::get_instance() {
        if (!instance) {
                instance = new IONameHandler();
                instance->initialize_root_dir();
                instance->initialize_placeholder();
        }
        return instance;
}
/**
 * TODO: Treat condition - error message?
 */
void IONameHandler::initialize_root_dir() {
        char cCurrentPath[FILENAME_MAX];

        if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
        {
            //return errno;
        }
        cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
        this->root_dir = cCurrentPath;
        this->initialize_output_dir();
}
void IONameHandler::initialize_output_dir() {
        PetscBool flg; //PetscBool instead of PetscTruth data type
        char dir[PETSC_MAX_PATH_LEN];     /* directory name */

        /* Initialize output directory path */
        PetscOptionsHasName(PETSC_NULL,"-o",&flg);
        if (flg) {
                PetscOptionsGetString(PETSC_NULL,"-o",dir,PETSC_MAX_PATH_LEN,&flg);
            this->output_dir = dir;
        } else {
                this->output_dir = this->root_dir;
        }
}

std::string IONameHandler::get_root_dir() {
        return this->root_dir;
}

std::string IONameHandler::get_output_dir() {
        return this->output_dir;
}

std::string IONameHandler::get_input_file_name(std::string file_name) {
        std::string file = this->get_root_dir() + "/" + file_name;
        return substitute_value(file);
}

std::string IONameHandler::get_output_file_name(std::string file_name) {
        std::string file = this->get_output_dir() + "/" + file_name;
        return substitute_value(file);
}

std::string IONameHandler::substitute_value(std::string file) {
        for (std::map<std::string,std::string>::const_iterator it = this->placeholder.begin(); it != this->placeholder.end(); ++it) {
                size_t i = file.find((*it).first,0);
                if(i != std::string::npos) {
                        file.replace(i, (*it).first.size(), (*it).second);
                }
        }
        return file;
}

bool IONameHandler::add_placeholder_item(std::string key, std::string value) {
        this->placeholder.insert( pair<std::string,std::string>(key,value));
        return true;
}

/*std::string IONameHandler::remove_placeholder_item(std::string key) {
        std::string ret = "";
        for (std::map<std::string,std::string>::iterator it = this->placeholder.begin(); it != this->placeholder.end(); it++) {
                if ((*it).first.compare(key) == 0) {
                        ret += (*it).second;
                        this->placeholder.erase(it);
                        break;
                }
        }
        return ret;
}
*/
void IONameHandler::initialize_placeholder() {
        PetscBool flg; // PetscBool instead of PetscTruth data type
        char dir[PETSC_MAX_PATH_LEN];     /* directory name */

        /* Initialize input directory name */
        PetscOptionsHasName(PETSC_NULL,"-i",&flg);
        if (flg) {
                PetscOptionsGetString(PETSC_NULL,"-i",dir,PETSC_MAX_PATH_LEN,&flg);
            this->add_placeholder_item("${INPUT}",dir);
        }
}
