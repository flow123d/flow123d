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
 * @file    read_ini.h
 * @brief   
 */

#ifndef READ_INI_H
#define READ_INI_H

#include <vector>

struct Ini_item;
struct Read_ini;

#include "system/system.hh"

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
	/**
	 * Read value of particular key as list of doubles and store them  into array.
	 * defval string is used if the key is not found.
	 */
	void	OptGetDblArray(const char *section, const char *key, const char *defval, std::vector<double> &array);

	void	OptGetIntArray(const char *section, const char *key, const char *defval, int Arrsize, int *Array);
	//char * OptGetStrArray();

#endif

