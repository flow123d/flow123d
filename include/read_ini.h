#ifndef READ_INI_H
#define READ_INI_H

struct Ini_item;
struct Read_ini;

#define DIR_DELIMITER '/'

struct Read_ini
{
	char *ini_file;
	char *ini_dir;     /* absolute or relative path to the ini file */
	struct Ini_item *ini_item;
};

struct Ini_item
{
	struct Ini_item *next;
	struct Ini_item *prev;
	char *section;
	char *key;
	char *value;

};
	long int OptGetInt(const char *section,const char *key,const char *defval);
	char *   OptGetStr(const char *section,const char *key,const char *defval);
	char *   OptGetFileName(const char *section,const char *key,const char *defval);
	bool     OptGetBool(const char *section,const char *key,const char *defval);
	double   OptGetDbl(const char *section,const char *key,const char *defval);
	void     OptionsInit(const char *fname );

#endif

