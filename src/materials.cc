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
 * @brief Read and initialize materials.
 *
 */

#include "system/system.hh"
#include <boost/algorithm/string/trim.hpp> // trim string leading and ending spaces
#include <boost/algorithm/string/classification.hpp>

//#include "mesh/mesh.h"
//#include "problem.h"
#include "materials.hh"
//#include "transport.h"
#include "system/math_fce.h"

flow::VectorId<int> xx;

/**
 * Default material constructor, set default values.
 */
Material::Material()
: dimension(0),
  hydrodynamic_resistence(NULL),
  por_m(1.0),
  por_imm(0.0),
  alpha(NULL),
  sorp_type(NULL),
  sorp_coef(NULL),
  size(1.0),
  phi(1.0),
  density(1.0),
  stor(1.0)
{}

/**
 * Default material destructor.
 */
Material::~Material()
{
    if (NONULL(hydrodynamic_resistence)) delete[] hydrodynamic_resistence;
    if (NONULL(alpha)) delete[] alpha;
    if (NONULL(sorp_type)) delete[] sorp_type;
    if (NONULL(sorp_coef)) delete[] sorp_coef;
}


/*!
 * @brief READ DATA OF ALL MATERIALS
 * @param file_name
 */
MaterialDatabase::MaterialDatabase(const string &file_name)
: file_name(file_name)
{
	FILE   *fin;   // input file
	char   line[ LINE_SIZE ];   // line of data file
	Iter mat;
	bool found;

    F_ENTRY;

    lock=false;
    xprintf( Msg, "Reading materials... ")/*orig verb 2*/;
    fin = xfopen( file_name.c_str(), "rt" );

    // read section with hydraulic conductivity tensor
    // create all new materials in the database
    found=skip_to( fin, "$Materials");
    INPUT_CHECK( found , "Can not find section: Materials.\n" );

    xfgets( line, LINE_SIZE - 2, fin );

    int n_materials = atoi( xstrtok(line) );
    INPUT_CHECK( n_materials >= 1 ,"Number of materials < 1 in section Materials\n");
    this->reserve(n_materials);

    while( n_materials-- > 0 ) {
        xfgets( line, LINE_SIZE - 2, fin );
        parse_and_add_material_line( line );
    }

    valid_sections["$Materials"]=this->size();
    xprintf( MsgVerb, " %d materials readed. ", this->size() )/*orig verb 4*/;
    xprintf(Msg, "O.K.\n");

    lock=true;

    // read cross cut of fractures
    xrewind(fin);
    found = skip_to(fin, "$Geometry");
    std::string cut_white;
    if (found) {
        int count=0;
        while (1) {
            xfgets(line, LINE_SIZE - 2, fin);
            cut_white=line;
            boost::algorithm::trim_if(cut_white, boost::algorithm::is_any_of(" \t\r\n")); //cut white space
            if ( cut_white.compare("$EndGeometry") == 0) break;

            int mid      = atoi( xstrtok(line, 0) );
            int dim      = atoi( xstrtok(NULL, 1) );
            double value    = atof( xstrtok(NULL, 2) );

            mat = find_id(mid);
            INPUT_CHECK( mat != this->end(), "Unknown material #%d found in section Geometry\n", mid);
            INPUT_CHECK( mat->dimension == dim, "Material dimension mismatch\n");
            mat->size=value;
            count++;
        }
        valid_sections["Geometry"] = count;
    }
    xfclose( fin );

    // read storativity for unsteady saturated flow
    read_one_parameter_section("Storativity", offsetof(Material, stor));
}

void MaterialDatabase::read_transport_materials(bool dual_porosity, bool sorption, int n_substances) {
    FILE *fin; // input file
    char line[LINE_SIZE]; // line of data file
    Iter mat;
    bool found;
    int sorption_n_coef_for_types[] = { 1, 2, 2 };
    std::string cut_white;

    F_ENTRY;

    fin = xfopen(file_name.c_str(), "rt");

    //if (dual_porosity) {   - have to be read for all transport problems, not only for dual porosity problems
        // read dual porosity
        found = skip_to(fin, "$DualPorosity");
        if (found) {
            int count = 0;
            while (1) {
                xfgets(line, LINE_SIZE - 2, fin);
                cut_white = line;
                boost::algorithm::trim_if(cut_white, boost::algorithm::is_any_of(" \t\r\n")); //cut white space
                if (cut_white.compare("$EndDualPorosity") == 0)
                    break;

                int mid = atoi(xstrtok(line, 0));
                mat = find_id(mid);
                INPUT_CHECK( mat != this->end(), "Unknown material #%d find in section %s\n", mid, "DualPorosity");

                mat->por_m = atof(xstrtok(NULL));

                if (dual_porosity) {
                    mat->por_imm = atof(xstrtok(NULL));
                    mat->alpha = new double[n_substances];
                    for (int ci = 0; ci < n_substances; ci++)
                        mat->alpha[ci] = atof(xstrtok(NULL));
                }
                count++;
            }
            valid_sections["DualPorosity"] = count;
        }
   // }

    if ( sorption ) {
	    xrewind(fin);
        found = skip_to(fin, "$Sorption");
        if (found) {
            // allocate  space for sorption parameters for substances
            for (Iter i = begin(); i != end(); i++) {
                i->sorp_type = new int[n_substances];
                i->sorp_coef = new std::vector<double>[n_substances];
            }

            int count = 0;
            while (1) {
                xfgets(line, LINE_SIZE - 2, fin);
                cut_white = line;
                boost::algorithm::trim_if(cut_white, boost::algorithm::is_any_of(" \t\r\n")); //cut white space
                if (cut_white.compare("$EndSorption") == 0)
                    break;

                int mid = atoi(xstrtok(line, 0));
                mat = find_id(mid);
                INPUT_CHECK( mat != this->end(), "Unknown material #%d find in section %s\n", mid, "Sorption");

                int si = atoi(xstrtok(NULL, 1)); // read substance index
                int type = mat->sorp_type[si] = atoi(xstrtok(NULL, 2));
                INPUT_CHECK(type >= 1 && type <= 3, "Wrong sorption type material #%d substance #%d.",this->get_id(mat),si);
                int n_sorp_coef = sorption_n_coef_for_types[type - 1];

                mat->sorp_coef[si].resize(n_sorp_coef);
                for (int ci = 0; ci < n_sorp_coef; ci++)
                    mat->sorp_coef[si][ci] = atof(xstrtok(NULL));

                count++;

            }
            valid_sections["Sorption"] = count / n_substances;
        }
    }
    xfclose(fin);

    // read sorption fraction in the case of dual_porosity
    if (dual_porosity) {
        read_one_parameter_section("SorptionFraction", offsetof(Material, phi));
    }
}

/*!
 * @brief  CREATE NEW MATERIAL
 * @return material structure
 */
MaterialDatabase::FullIter MaterialDatabase::new_material( int id )
{

	F_ENTRY;

    if (lock) xprintf(PrgErr, "Try to add new material into locked database!");

	return( add_item(id));
}


// *****************************************************************************************


#define K_1D      11
#define A_1D     -11
#define K_2D_1    21
#define K_2D_2    22
#define K_2D_3    23
#define A_2D_1   -21
#define A_2D_2   -22
#define A_2D_3   -23
#define K_3D_1    31
#define K_3D_3    33
#define K_3D_6    36
#define A_3D_1   -31
#define A_3D_3   -33
#define A_3D_6   -36


/*!
 * @brief SET THE "a[][]" AND "k[][]" FIELDS IN STRUCT ELEMENT
 * @param[in,out] ele element
 */
void MaterialDatabase::calc_material_resistance( int type, double coef[6], MaterialDatabase::Iter mat)
{
    F_ENTRY;

    // allocate always full 3D tensor, initialize by zero
    mat->hydrodynamic_resistence = new double[9];

    SmallMtx1 mtx1=(SmallMtx1) mat->hydrodynamic_resistence;
    SmallMtx2 mtx2=(SmallMtx2) mat->hydrodynamic_resistence;
    SmallMtx3 mtx3=(SmallMtx3) mat->hydrodynamic_resistence;

    switch( type ) {
        case K_1D:
        case A_1D:
            mtx1[ 0 ][ 0 ] = coef[ 0 ];
            ASSERT(mat->dimension == 1, "Material dimension mismatch\n");
            break;
        case K_2D_1:
        case A_2D_1:
            coef[ 1 ] = coef[ 0 ];
        case K_2D_2:
        case A_2D_2:
        case K_2D_3:
        case A_2D_3:
            mtx2[ 0 ][ 0 ] = coef[ 0 ];
            mtx2[ 0 ][ 1 ] = coef[ 2 ];
            mtx2[ 1 ][ 0 ] = coef[ 2 ];
            mtx2[ 1 ][ 1 ] = coef[ 1 ];
            ASSERT(mat->dimension == 2, "Material dimension mismatch\n");
            break;
        case K_3D_1:
        case A_3D_1:
            coef[ 1 ] = coef[ 0 ];
            coef[ 2 ] = coef[ 0 ];
        case K_3D_3:
        case A_3D_3:
        case K_3D_6:
        case A_3D_6:
            mtx3[ 0 ][ 0 ] = coef[ 0 ];
            mtx3[ 0 ][ 1 ] = coef[ 3 ];
            mtx3[ 0 ][ 2 ] = coef[ 4 ];
            mtx3[ 1 ][ 0 ] = coef[ 3 ];
            mtx3[ 1 ][ 1 ] = coef[ 1 ];
            mtx3[ 1 ][ 2 ] = coef[ 5 ];
            mtx3[ 2 ][ 0 ] = coef[ 4 ];
            mtx3[ 2 ][ 1 ] = coef[ 5 ];
            mtx3[ 2 ][ 2 ] = coef[ 2 ];
            ASSERT(mat->dimension == 3, "Material dimension mismatch\n");
            break;
    }

    if( type > 0 ) {
        // if conductivity is given, we have to calculate the inverse
        double *conductivity=mat->hydrodynamic_resistence;
        mat->hydrodynamic_resistence = new double[9];

        double det= MatrixInverse(conductivity,mat->hydrodynamic_resistence, mat->dimension);
        if ( DBL_EQ( 0.0, det ) )
            xprintf(UsrErr,"Conductivity of material id: %d is almost singular. Determinant: %g \n", get_id(mat),det);

        delete[] conductivity;
    }
}


/*!
 * @brief PARSE hydraulic conductivity / resistence MATERIAL LINE
 * @param[in]  line input text line
 */
void MaterialDatabase::parse_and_add_material_line(char *line )
{
	F_ENTRY;

	ASSERT(NONULL(line),"NULL as argument of function parse_material_line()\n");

	int mid   = atoi( xstrtok(line,0) );
    INPUT_CHECK( mid >= 0 ,"Id number of material must be at least 0.\n");
	Iter mat = new_material(mid);

	int type  = atoi( xstrtok(NULL,1) );
	// check correct type, get dimension
    switch( type ) {
	        case K_1D:
	        case A_1D:
                mat->dimension=1;
                break;
	        case K_2D_1:
	        case K_2D_2:
	        case K_2D_3:
	        case A_2D_1:
	        case A_2D_2:
	        case A_2D_3:
                mat->dimension=2;
                break;
            case K_3D_1:
	        case K_3D_3:
	        case K_3D_6:
	        case A_3D_1:
	        case A_3D_3:
	        case A_3D_6:
                mat->dimension=3;
                break;
	        default:
	            INPUT_CHECK( false , "Wrong resistance input type %d\n", type);
	}
    // read coeficiants
    double coef[6];
	SET_ARRAY_ZERO(coef,6);

	int n_coef=abs(type) % 10;
	for(unsigned int  ci = 0; ci < n_coef; ci++ )
		coef[ ci ]  = atof( xstrtok(NULL) );

	calc_material_resistance(type, coef, mat);
}






/**
 * read a section given by @param section_name without '$'
 * fill values into struct Material into a double variable with offset given by @param offset_to_set
 */
void MaterialDatabase::read_one_parameter_section(const string &section_name, const uint8_t offset_to_set)
{
    F_ENTRY_P(section_name);

    FILE   *fin;   // input file
    char   line[ LINE_SIZE ],string[ LINE_SIZE ];   // line of data file
    bool found;

    std::string start_str("$"); start_str+=section_name;
    std::string end_str("$End"); end_str+=section_name;
    std::string cut_white;

    fin = xfopen( file_name.c_str(), "rt" );
    found = skip_to(fin, start_str.c_str());
    if (found) {
        int count=0;
        while (1) {
            xfgets(line, LINE_SIZE - 2, fin);
            cut_white=line;
            boost::algorithm::trim_if(cut_white, boost::algorithm::is_any_of(" \t\r\n")); //cut white space
            if ( cut_white.compare(end_str) == 0) break;

            int mid      = atoi( xstrtok(line, 0) );
            double value    = atof( xstrtok(NULL, 1) );

            Iter mat = find_id(mid);
            INPUT_CHECK( mat != this->end(), "Unknown material #%d found in section %s\n", mid, section_name.c_str());
            double *entry_to_set =(double *)( (uint8_t *)(&(*mat)) + offset_to_set );
            *entry_to_set = value;
            count++;
        }
        valid_sections[section_name] = count;
    }
    xfclose( fin );
}



//-----------------------------------------------------------------------------
// vim: set cindent:









