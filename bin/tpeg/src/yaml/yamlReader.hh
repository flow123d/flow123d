
#include <iostream>

#include <yaml.h>
#include <vector>

using namespace std;

class YamlReader {
   private:
      const char *filename;
      FILE *ymlSettings;
      yaml_parser_t parser;
      yaml_event_t input_event;
      int walkLevel(vector<string> xPath, vector<string> *result, int outputEvent, bool res = false);
      int outputEvent(int inputEvent);
      int initialize();

   public:
      YamlReader(const char * __filename);
      ~YamlReader();

      const char * getFilename();
      int getValue(vector<string> xPath, vector<string> *result);
      int getValue(vector<string> xPath, string *result);
      int getValue(vector<string> xPath, vector<double> *result);
      int getValue(vector<string> xPath, double *result);
      int getValue(vector<string> xPath, vector<int> *result);
      int getValue(vector<string> xPath, int *result);
};
