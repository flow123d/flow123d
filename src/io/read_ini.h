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
 * @file
 * @brief ???
 *
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

