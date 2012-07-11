/*
 * readTest.cc
 *
 *  Created on: Nov 23, 2010
 *      Author: honzaL
 */

#include<iostream>
#include<yaml.h>
#include<vector>

#include"yamlReader.hh"

using namespace std;

YamlReader::YamlReader(const char * __filename) {
    filename = __filename;

    /* Opening init file */

    ymlSettings = fopen(filename, "rb");
    if (ymlSettings == NULL) {
        fprintf(stderr, "Error opening file %s", filename);
        exit(1);
    }

    /* Clear the object. */

    memset(&parser, 0, sizeof (parser));

    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "Could not initialize the parser object\n");
        exit(1);
    }

    /* Set the parser parameters. */

    yaml_parser_set_input_file(&parser, ymlSettings);
}

YamlReader::~YamlReader() {
    if (ymlSettings != NULL) {
        fclose(ymlSettings);
    }
    if (&parser != NULL) {
        yaml_parser_delete(&parser);
    }
}

int YamlReader::initialize() {
    if (&parser != NULL) {
        yaml_parser_delete(&parser);
    }

    /* Clear the object. */

    memset(&parser, 0, sizeof (parser));

    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "Could not initialize the parser object\n");
        exit(1);
    }

    /* Set the parser parameters. */
    rewind(ymlSettings);
    yaml_parser_set_input_file(&parser, ymlSettings);

    if (!yaml_parser_parse(&parser, &input_event)) {
        fprintf(stderr, "EXIT in YamlReader::initialize();\n");
        exit(1);
    }
    return 1;
}

int YamlReader::getValue(vector<string> xPath, int *result) {
    initialize();
    vector<string> newPath;
    int i = 0;
    for (vector<string>::iterator xPathIterator = (xPath.end() - 1); xPathIterator != xPath.begin() - 1; xPathIterator--, i++) {
        newPath.resize(i + 1, *xPathIterator);
    }
    vector<string> res;
    int ret = walkLevel(newPath, &res, outputEvent(input_event.type));
    if (ret) {
        (*result) = atoi(res[0].c_str());
    }
    return ret;
}

int YamlReader::getValue(vector<string> xPath, vector<int> *result) {
    initialize();
    vector<string> newPath;
    int i = 0;
    for (vector<string>::iterator xPathIterator = (xPath.end() - 1); xPathIterator != xPath.begin() - 1; xPathIterator--, i++) {
        newPath.resize(i + 1, *xPathIterator);
    }
    vector<string> res;
    int ret = walkLevel(newPath, &res, outputEvent(input_event.type));
    if (ret) {
        result->resize(res.size(), 0);
        for (int i = 0; i < res.size(); i++) {
            (*result)[i] = atoi(res[i].c_str());
        }
    }
    return ret;
}

int YamlReader::getValue(vector<string> xPath, double *result) {
    initialize();
    vector<string> newPath;
    int i = 0;
    for (vector<string>::iterator xPathIterator = (xPath.end() - 1); xPathIterator != xPath.begin() - 1; xPathIterator--, i++) {
        newPath.resize(i + 1, *xPathIterator);
    }
    vector<string> res;
    int ret = walkLevel(newPath, &res, outputEvent(input_event.type));
    if (ret) {
        (*result) = atof(res[0].c_str());
    }
    return ret;
}

int YamlReader::getValue(vector<string> xPath, vector<double> *result) {
    initialize();
    vector<string> newPath;
    int i = 0;
    for (vector<string>::iterator xPathIterator = (xPath.end() - 1); xPathIterator != xPath.begin() - 1; xPathIterator--, i++) {
        newPath.resize(i + 1, *xPathIterator);
    }
    vector<string> res;
    int ret = walkLevel(newPath, &res, outputEvent(input_event.type));
    if (ret) {
        result->resize(res.size(), 0);
        for (int i = 0; i < res.size(); i++) {
            (*result)[i] = atof(res[i].c_str());
        }
    }
    return ret;
}

/**
 * 
 * @param xPath
 * @param result
 * @return 
 */
int YamlReader::getValue(vector<string> xPath, string *result) {
    //cout << "      - YamlReader::getValue(vector<string>, string*) - Start" << endl;

    //cout << "         - Start initialization" << endl;
    initialize();
    //cout << "         - End initialization" << endl;

    //cout << "         - Start newPath" << endl;
    vector<string> newPath;
    int i = 0;
    for (vector<string>::iterator xPathIterator = (xPath.end() - 1); xPathIterator != xPath.begin() - 1; xPathIterator--, i++) {
        newPath.resize(i + 1, *xPathIterator);
    }
    //cout << "         - End newPath" << endl;

    //cout << "         - Start result" << endl;
    vector<string> res;
    int ret = walkLevel(newPath, &res, outputEvent(input_event.type));
    if (ret) {
        (*result) = res[0];
    }
    //cout << "         - End result" << endl;

    //cout << "      - YamlReader::getValue(vector<string>, string*) - End" << endl;

    return ret;
}

int YamlReader::getValue(vector<string> xPath, vector<string> *result) {
    initialize();
    //	vector<string> xPath(xPathArr, xPathArr + sizeof(xPathArr) / sizeof(string));
    vector<string> newPath;
    int i = 0;
    for (vector<string>::iterator xPathIterator = (xPath.end() - 1); xPathIterator
            != xPath.begin() - 1; xPathIterator--, i++) {
        newPath.resize(i + 1, *xPathIterator);
    }
    return walkLevel(newPath, result, outputEvent(input_event.type));
}

const char * YamlReader::getFilename() {
    return filename;
}

int YamlReader::outputEvent(int inputEvent) {
    int ret;
    switch (inputEvent) {
        case YAML_DOCUMENT_START_EVENT:
            ret = YAML_DOCUMENT_END_EVENT;
            break;
        case YAML_STREAM_START_EVENT:
            ret = YAML_STREAM_END_EVENT;
            break;
        case YAML_SEQUENCE_START_EVENT:
            ret = YAML_SEQUENCE_END_EVENT;
            break;
        case YAML_MAPPING_START_EVENT:
            ret = YAML_MAPPING_END_EVENT;
            break;
        default:
            ret = YAML_NO_EVENT;
    }
    if (ret == YAML_NO_EVENT) {
        cout << "EXIT in YamlReader::outputEvent(int inputEvent);" << endl;
        exit(1);
    }
    return ret;
}

int YamlReader::walkLevel(vector<string> xPath, vector<string> *result, int outputEvent, bool res) {
    int done = 0;
    int ret = 0;
    bool tmp = false;
    vector<string> newPath;
    bool key = (outputEvent == YAML_SEQUENCE_END_EVENT) ? false : true;
    while (!done) {

        if (!yaml_parser_parse(&parser, &input_event)) {
            cout << "EXIT in vector<string> xPath, vector<string> *result, int outputEvent, bool res);" << endl;
            exit(1);
        }
        done = (input_event.type == outputEvent);

        // Analyze the event.

        switch (input_event.type) {
            case YAML_MAPPING_START_EVENT:
            case YAML_SEQUENCE_START_EVENT:
            case YAML_DOCUMENT_START_EVENT:
            case YAML_STREAM_START_EVENT:
                if (outputEvent == YAML_DOCUMENT_END_EVENT || outputEvent == YAML_STREAM_END_EVENT) {
                    done += ret += walkLevel(xPath, result, this->outputEvent(input_event.type), res);
                } else {
                    done += ret += walkLevel(newPath, result, this->outputEvent(input_event.type), res);
                }
                res = false;
                key = true;
                break;
            case YAML_SCALAR_EVENT:
                if (key) {
                    if (!xPath.empty() && xPath.back() == (char *) input_event.data.scalar.value) {
                        // jsem ve spravne vetvi
                        if (xPath.size() > 1) { // muzu kopirovat cestu, nejsem v koncovem nodu
                            newPath.resize(xPath.size() - 1);
                            copy(xPath.begin(), xPath.end() - 1, newPath.begin());
                        } else if (xPath.size() == 1) { // dalsi hodnota (scalar nebo sequence of scalar) je hledana hodnota
                            res = true;
                        }
                        key = false;
                    } else {
                        res = false;
                    }
                } else { // not key
                    if (res) {
                        result->resize(result->size() + 1, (char *) input_event.data.scalar.value);
                        ret = 1;
                    }
                    if (res && outputEvent == YAML_MAPPING_END_EVENT) { // zpatecni je mapovaci, tedy jde pouze o jednu hodnotu
                        res = false;
                    }

                }
                break;

            default:
                break;
        }
    }
    return ret;
}
