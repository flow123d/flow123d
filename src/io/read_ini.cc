/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    read_ini.cc
 * @ingroup io
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

#include "system/system.hh"
#include "system/xio.h"
#include "io/read_ini.h"

#include <boost/tokenizer.hpp>
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>

static struct Read_ini *read_ini = NULL;

#define FOR_INI_ITEMS(i)     for((i)=read_ini->ini_item;(i)!=NULL;(i)=(i)->next)

static void make_ini_item_list(const char *fname);
static char *section_test(char *section);
static char *strip_spaces(char *string);
static struct Ini_item *new_item(struct Ini_item *prev,char *section, char *key, char *value);

//#define	xOptGetStr(s,k)	(ini.GetValue(s,k,NULL))




/*!
 * @brief      STRTOK WITH ERROR HANDLING and whitespace delimiters
 * @param[in]  s        strtok string pointer
 * @param[in]  position requested position of the token
 * @return              strtok return
 */
char *xstrtok(char *s, int position)
{
    char *rc;
    const char * const whitespace_delim=" \t\r\n";

    rc = xstrtok( s, whitespace_delim, position);
    return(rc);
}

/*!
 * @brief      STRTOK WITH ERROR HANDLING and user specified delimiters
 * @param[in]  s1       strtok string pointer
 * @param[in]  delim    delimiters
 * @param[in]  position requested position of the token
 * @return              strtok return
 *
 * Function behaves like original strtok
 */
char *xstrtok( char *s1, const char *delim, int position )
{
    char *rc;
    static char * full_string = NULL;
    static int token_count;

    OLD_ASSERT(!( delim == NULL ),"NULL pointer as delimiter in xstrtok()\n");

    if ( s1 )
    {
        if ( !full_string )
        {
            full_string = (char *)xmalloc( LINE_SIZE );
            full_string[0] = 0x0;
        }

        strncpy( full_string, s1, LINE_SIZE );
        token_count = 0;
    }

    INPUT_CHECK( token_count == position || position < 0, "Requested position %d dosn't match token position %d", position, token_count);
    rc = strtok( s1, delim );
    token_count++;

    INPUT_CHECK(!( rc == NULL ),"Missing token no. %d: original string '%s' with delimiters '%s'\n", token_count, full_string, delim );

    return(rc);
}



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

	read_ini=(struct Read_ini*)xmalloc(sizeof(struct Read_ini));

	ini=xfopen(fname,"rt");
	ASSERT(ini)(fname).error("Failed to open the ini file");


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
    xprintf(Err, "OptGetXXX input interface is not supported anymore.\n");


	const char *rc = NULL;
	struct Ini_item *ini_item;

	FOR_INI_ITEMS(ini_item)
		if( (!strcmp(ini_item->section,section)) && (!strcmp(ini_item->key,key)) ){
			rc = ini_item->value;
			break;
		}

	if (rc == NULL) {
		if (defval == NULL)
			xprintf(UsrErr,"Required parameter: section '%s' key '%s' is not given.\n",section,key);
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
bool OptGetBool( const char *section,const  char *key,const  char *defval )
{
	char *str;
	char res=false;
	str = OptGetStr(section, key, defval);

	if ( boost::iequals(str, "yes") || boost::iequals(str, "true") || boost::iequals(str, "1") ) res=true;
	else if ( boost::iequals(str, "no") || boost::iequals(str, "false") || boost::iequals(str, "0") ) res=false;
	else {
		xfree(str);
		if (defval == NULL) xprintf(UsrErr,"Required parameter: [%s] %s is not a boolen.\n",section,key);
		str=(char *)defval;
		if ( boost::iequals(str, "yes") || boost::iequals(str, "true") || boost::iequals(str, "1") ) res=true;
		else if ( boost::iequals(str, "no") || boost::iequals(str, "false") || boost::iequals(str, "0") ) res=false;
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

	ASSERT(fname).error("NULL file name.\n");

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

//=============================================================================
// GET DOUBLE ARRAY VARIABLE FROM INI FILE
//=============================================================================
void OptGetDblArray( const char *section,const  char *key,const  char *defval, std::vector<double> &array)
// ArrSize contain number of Array members (length), *Array is the adress of array which should be filled up
{

	char * tmp_str = OptGetStr(section,key,defval);
	std::string str = tmp_str;
	free(tmp_str);
	boost::tokenizer<boost::char_separator<char> > line_tokenizer(str, boost::char_separator<char>("\t "));
	boost::tokenizer<boost::char_separator<char> >::iterator tok;

	double value;
	try {
	    for(    tok = line_tokenizer.begin();
	            tok != line_tokenizer.end();
	            ++tok) {
	        value = boost::lexical_cast<double> (*tok);
	        array.push_back(value);
	    }
	}  catch (boost::bad_lexical_cast &) {
        xprintf(UsrErr, "INI file: Can not convert token `%s` of key `[%s] %s` to double.\n", (*tok).c_str(), section, key);
    }
}


//=============================================================================
// GET Int ARRAY VARIABLE FROM INI FILE
//=============================================================================
void OptGetIntArray( const char *section,const  char *key,const  char *defval, int ArrSize, int *Array)// ArrSize contain number of Array members (length), *Array is the adress of array which should be filled up
{
	char *str;
	int res;
	int i;

	str=OptGetStr(section,key,defval);
	for(i = 1; i < ArrSize; i++){
		if (sscanf(str,"%d",&res) == 0) {
		if (defval == NULL) xprintf(UsrErr,"Can not convert %d. ini-file entry to integer parameter: [%s] %s.\n",i,section,key);
		if (sscanf(defval,"%d",&res) == 0)
			xprintf(PrgErr,"Default value \"%s\" of parameter: [%s] %s is not an integer.\n",defval,section,key);
		}else{
		  *(Array + (i-1)*sizeof(double)) = res;
		}
	}

	free( str );
	return;
}

//=============================================================================
// GET string ARRAY VARIABLE FROM INI FILE
//=============================================================================
/*char *OptGetStrArray(const char *section,const char *key, int sb_count, struct TS_lat *dest)
{
	const char **rc = NULL;
	struct Ini_item *ini_item;
	int i;

	FOR_INI_ITEMS(ini_item)
		if( (!strcmp(ini_item->section,section)) && (!strcmp(ini_item->key,key)) ){
			for(i=0; i < sb_count; i++)
			{
			  if(sscanf(ini_item->value,"%s",dest[i].nazev) == NULL)
			  {
				printf("\nerror during required %d-th parameter initialization occured\n",i);
			  }else{
			  	printf("\nthe name of %d-th substance is %s\n",i,dest[i].nazev);
			  }
			}
			//
			*rc = ini_item->value;
			break;
		}

	if (rc == NULL) {
		if (defval == NULL)
			xprintf(UsrErr,"Required parameter: [%s] %s is not given.\n",section,key);
		else
			*rc = defval;
	}
	return xstrcpy(*rc);
}*/

