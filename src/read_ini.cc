/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 *
 * @file
 * @brief   OPTIONS RUTINES - get program parameters, reading from options/ini file
 * @section DESCRIPTION
 *
 *   Functions OptGet* are implemented through OptGetStr - only one interface routine for
 *   various implementation the parameter default is always char *, which can be NULL in the case that
 *   explicit value is necessary and program should produce an error if corresponding
 *   parameter is not given.
 *
 *   Current implementation is based on SimpleINI C++ library.
 *
 */

#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

#include "system.hh"
#include "xio.h"
#include "read_ini.h"

static struct Read_ini *read_ini = NULL;

#define FOR_INI_ITEMS(i)     for((i)=read_ini->ini_item;(i)!=NULL;(i)=(i)->next)

static void make_ini_item_list(const char *fname);
static char *section_test(char *section);
static char *strip_spaces(char *string);
static struct Ini_item *new_item(struct Ini_item *prev,char *section, char *key, char *value);

//#define	xOptGetStr(s,k)	(ini.GetValue(s,k,NULL))

//=============================================================================
// MAKE INI KEYS LIST
//=============================================================================
void make_ini_item_list(const char *fname)
{
	FILE *ini;

	char line[ LINE_SIZE ];
	char string[ LINE_SIZE ];
	char section[ LINE_SIZE ];//="";
	char *section_ptr;
	char *key;
	char *value;
	char *tmp;
	struct Ini_item *prev = NULL;
	int i;

	read_ini=(struct Read_ini*)xmalloc(sizeof(struct Read_ini));

	ini=xfopen(fname,"rt");
	ASSERT(NONULL(ini),"Failed to open the ini file: %s",fname);


	while( xfgets( line, LINE_SIZE - 2, ini ) != NULL ) {
	    sscanf( line, "%s", string );      // store first substring in the string


	    if (strlen(string)==0)             // skips start blank lines
	        continue;

		// READ SECTION
	    section_ptr = section_test(string);
		if ( section_ptr ) { // test of new section
		    strcpy(section,section_ptr);
			continue;					// go to next line
		}

		if(strchr(line,'=') == NULL)	// test of "=" character on the line
			continue;

		// CLEAR COMMENTS
		tmp = strstr(line,"#"); 	// find comment position
		if(tmp != NULL){
			sprintf(tmp,"%s","");	// force ends string on comment position
			if(strlen(line) == 0)	// continue if line contains only comment
				continue;			// go to next line
		}

		// READ KEY
		tmp = xstrtok(line,"=");	// read characters before "="
		sscanf(tmp,"%s",string);	// strip spaces
		if(strlen(string) == 0)
			continue;
		else
			key = xstrcpy(string);


		//READ VALUE
		tmp = xstrtok(NULL,"=");
		tmp = strip_spaces(tmp);
		if(strlen(tmp) == 0){ //string
			xfree(key);
			continue;
		}
		else
			value = xstrcpy(tmp);

/*
		printf("%s\n",section);
		printf("%s\n",key);
		printf("%s\n\n",value);
*/

		prev = new_item(prev,section,key,value);
		xfree(key);
		xfree(value);
	};
};
//=============================================================================
// STRIP START AND END BLANK CHARACTERS
//=============================================================================
char *strip_spaces(char *string)
{
	int i;
	while((string[0] ==' ') || (string[0] =='\t')){
		string++;
	}
	i = strlen(string) - 1;
	while((string[i] ==' ') || (string[i] =='\t') || (string[i] =='\n')  || (string[i] =='\r')){
		string[i--] = 0;
	}
	return string;
}
//=============================================================================
// ADD NEW KEY TO INI-LIST
//=============================================================================
struct Ini_item *new_item(struct Ini_item *prev,char *section, char *key, char *value)
{
	struct Ini_item *ini_item;

	if((section != NULL) && (key != NULL) && (value != NULL)){
		ini_item=(struct Ini_item*)xmalloc(sizeof(struct Ini_item));

		if(prev == NULL){
			read_ini->ini_item = ini_item;
			ini_item->prev = NULL;
		}
		else{
			ini_item->prev = prev;
			ini_item->next = NULL;
			prev->next = ini_item;
		}

		ini_item->section = xstrcpy(section);
		ini_item->key = xstrcpy(key);
		ini_item->value = xstrcpy(value);
		return ini_item;
	}
	else
		return prev;

}
//=============================================================================
// TEST OF SECTION STRING
//=============================================================================
char *section_test(char *line)
{
	if(line == NULL) return NULL;

	if( (line[0] == '[') && (line[strlen(line)-1] == ']') && (strlen(line) > 2))
		return (xstrtok(line,"[]"));
	else
		return (NULL);
};

/*!
 * @brief Create new string from selected variable from ini file
 * @param[in] section  Inifile section
 * @param[in] key      Inifile key
 * @param[in] defval   Default value, if key is not found in inifile
 * @return Ptr to new string with option value from Inifile or to default value
 */
char *OptGetStr(const char *section,const char *key,const char *defval)
{
	const char *rc = NULL;
	struct Ini_item *ini_item;

	FOR_INI_ITEMS(ini_item)
		if( (!strcmp(ini_item->section,section)) && (!strcmp(ini_item->key,key)) ){
			rc = ini_item->value;
			break;
		}

	if (rc == NULL) {
		if (defval == NULL)
			xprintf(UsrErr,"Required parameter: [%s] %s is not given.\n",section,key);
		else
			rc = defval;
	}

	return xstrcpy(rc);
}

//=============================================================================
// GET FILE NAME VARIABLE FROM INI FILE
//=============================================================================
char *OptGetFileName(const char *section,const char *key,const char *defval)
{


    return OptGetStr(section,key,defval);
}

//=============================================================================
// GET INT VARIABLE FROM INI FILE
//=============================================================================
long int OptGetInt( const char *section,const char *key,const char *defval )
{
	char *str;
	long int res;

	str=OptGetStr(section,key,defval);
	if (sscanf(str,"%ld",&res) == 0) {
		if (defval == NULL) xprintf(UsrErr,"Can not convert to integer parameter: [%s] %s.\n",section,key);
		xprintf(PrgErr,"Default value %s of parameter: [%s] %s is not an integer.\n",defval,section,key);
	}

	xfree( str );
	return res;
}

//=============================================================================
// GET DOUBLE VARIABLE FROM INI FILE
//=============================================================================
double OptGetDbl( const char *section,const  char *key,const  char *defval )
{
	char *str;
	double res;

	str=OptGetStr(section,key,defval);
	if (sscanf(str,"%lg",&res) == 0) {
		if (defval == NULL) xprintf(UsrErr,"Can not convert to double parameter: [%s] %s.\n",section,key);
		if (sscanf(defval,"%lg",&res) == 0)
			xprintf(PrgErr,"Default value \"%s\" of parameter: [%s] %s is not an double.\n",defval,section,key);
	}

	xfree( str );
	return res;
}

//=============================================================================
// GET BOOL VARIABLE FROM INI FILE
//=============================================================================
#define SCMP(x)	(strcmpi( str, x ) == 0)
bool OptGetBool( const char *section,const  char *key,const  char *defval )
{
	char *str;
	char res;
	str = OptGetStr(section, key, defval);

	if ( SCMP("yes")||SCMP("true")||SCMP("1") ) res=true;
	else if ( SCMP("no")||SCMP("false")||SCMP("0") ) res=false;
	else {
		xfree(str);
		if (defval == NULL) xprintf(UsrErr,"Required parameter: [%s] %s is not a boolen.\n",section,key);
		str=(char *)defval;
		if ( SCMP("yes")||SCMP("true")||SCMP("1") ) res=true;
		else if ( SCMP("no")||SCMP("false")||SCMP("0") ) res=false;
		else xprintf(PrgErr,"Default value \"%s\" of parameter: [%s] %s is not a boolean.\n",defval,section,key);
	}
	return res;
}
/*!
 * @brief 			Load options file
 * @param[in] fname File name
 */
void OptionsInit(const char *fname )
{
	//char *path;
	//int len;

	ASSERT(NONULL(fname),"NULL file name\n");

	// take absolute path to the file
	// this is completly wrong in the case the absolute path is alredy given
	// since we should remain in the workdir when calling this it is not critical 
	// to hava correct function for path manipulation (see BOOST)
	/*
	path=xgetcwd();
	len=strlen(path)+strlen(fname)+1;
	if ( len > PATH_MAX )
	{
	    xprintf(UsrErr, "Path too long\n");
	}
	strcpy( options_fname, path ); xfree(path);
	strcat( options_fname, PATH_SEP );
	strcat( options_fname, fname );
*/

	// initialization of the INI reader
	make_ini_item_list(fname);
	
}
