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
 * $Id: test.cpp 1055 2011-04-21 13:43:54Z jan.brezina $
 * $Revision: 1055 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2011-04-21 15:43:54 +0200 (Thu, 21 Apr 2011) $
 *
 * @file test.cpp
 * @brief  Test case
 *
 */
/*
#include "system/system.hh"

int main(int argc, char * argv[]) {
  
  int *i = new int (2);
  INPUT_CHECK( *i == 2 , "new operator ... don't pass\n");
  delete i;
  
  int *ii = new int [10];
  INPUT_CHECK (ii[0] == 0, "new[] operator ... don't pass\n");
  delete ii;
  
  int (*iii)[3];
  iii=new (int[3]);
}  
*/