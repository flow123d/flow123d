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
 * @file   constantdb.h
 * @author dalibor
 * @date October 1, 2010, 10:55 PM
 */

#include <map>
#include <string.h>

#ifndef _CONSTANTDB_H
#define	_CONSTANTDB_H

#ifdef	__cplusplus
extern "C" {
#endif

    /*!
     *  @brief <B>DON'T USE THIS CLASS - it's just for my private using.</B>
     *
     * <H1>! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !</H1>
     * <H1>! !    DON'T USE THIS CLASS - it's just for my private using.     ! !</H1>
     * <H1>! !    After finishing my work,                                   ! !</H1>
     * <H1>! !     this class will be canceled without refund.               ! !</H1>
     * <H1>! !                                                      Dalibor  ! !</H1>
     * <H1>! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !</H1>
     * <H1></H1>
     * <H1>Design Pattern - Singleton</H1>
     * <BR>
     */
    class ConstantDB {
    private:
        static ConstantDB* instance;

        struct cmp_str {

            bool operator()(char const* a, char const* b) {
                return strcmp(a, b) < 0;
            }
        };

        std::map<const char*, const char*, cmp_str> charMap;
        std::map<const char*, int, cmp_str> intMap;
        std::map<const char*, double, cmp_str> doubleMap;
        std::map<const char*, void*, cmp_str> objectMap;

        /*!
         * @brief Deny construktor - instance is accesible by method getInstance().
         */
        ConstantDB() {
        }

        /*!
         * @brief Deny destruktor - instance died with the application.
         */
        ~ConstantDB() {
        }

        bool isKeyChar(const char*);
        bool isKeyInt(const char*);
        bool isKeyDouble(const char*);
        bool isKeyObject(const char*);

    public:

        /*!
         * @brief get instance (one and only one) of ConstantDB
         */
        static ConstantDB* getInstance() {
            return instance;
        }

        void setChar(const char*, const char*);
        const char* getChar(const char*);

        void setInt(const char*, int);
        int getInt(const char*);

        void setDouble(const char*, double);
        double getDouble(const char*);

        void setObject(const char*, void*);
        void* getObject(const char*);
    };


#ifdef	__cplusplus
}
#endif

#endif	/* _CONSTANTDB_H */

