/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Read and initialize materials.
 *
 * TODO: There should not be one structure gathering various physical values for various models.
 * Every model should have its own structures for its material data possibly with some linkage into a
 * global library of materials.
 *
 */
#include <map>
#include <vector>
#include <string>
#include <sys_vector.hh>


#ifndef MAKE_MATERIALS_H
#define MAKE_MATERIALS_H

/// @brief STRUCTURE OF THE MATERIAL, data read from material file
class Material
{

public:
    Material();
    ~Material();

//    int       id;           ///< input Id # of the material
//    int       order;        ///< internal continuous Id
    int       dimension;    ///< dimension of material (for tensor and vector values)
    double   *hydrodynamic_resistence;         ///< Values of coefficients
    double    por_m;        ///< Mobile porosity
    double    por_imm;      ///< Immobile porosity
    double   *alpha;        ///< Coefficient of nonequillibrum exchange
    int      *sorp_type;    ///< Type of sorption
    std::vector<double> *sorp_coef;    ///< Coefficient of sorption
    double    size;         ///< cross cut for materials of dimmension < 3 (for dim == 3 this is always == 1)
    double    phi;          ///< solid / solid mobile
    double    density;      ///< density of solid rock
    double    stor;         ///< Storativity of material
};


/**
 * This class provides access to the Material structure for given material Id. It also provides
 *  methods for reading meterial input file.
 *
 *  Current state is far to be optimal. Issues:
 *  1) More general reading methods
 *  2) Access to the values through run-time names not through structure.
 *  3) Possibility to create fast element-wise object for access to one scaler/vector/tensor value
 *  4) Dimensions.
 *  5) More general of special keywords on input.
 */

class MaterialDatabase : public flow::VectorId< Material > {

public:
    typedef flow::VectorId< Material >::Iter MaterialIter;

    /// Constructor. Reads a given file and fill the database with non-transport data.
    MaterialDatabase(const string &file_name);
    /**
     * Reads transport related sections. Unfortunately actual convention is that unimportant sections
     * are not read so they can contain parsing errors. In order to keep this behavior we have this
     * separate reading method. Moreover it needs number of transported substances.
     */
    void read_transport_materials(bool dual_porosity, bool sorption, int n_substances);

    /**
     * Lock the database after reading.
     */
    inline void lock_base(bool new_lock)
    {lock=new_lock;}

    MaterialIter new_material(int id);

    /// Returns true if the given section was read and have values for all materials.
    inline bool valid_section(const string &sec) {
        std::map<string, unsigned int>::iterator iter = valid_sections.find(sec);
        return ( iter != valid_sections.end() &&  iter->second == this->size() ); // check existence of section and its size
    }


    virtual ~MaterialDatabase() {}

private:

    /// material file name
    std::string file_name;
    /// std::map to store successfully reed sections and their sizes
    std::map<string,unsigned int> valid_sections;
    /// Indicates that no material should be added to the base.
    bool lock;

    /// Allocates and initialize new Material structure.
    void parse_and_add_material_line(char *line);

    /// Reads a section which sets just one scalar parameter.
    void read_one_parameter_section(const string &section_name, const uint8_t offest_to_set);

    /// Compute hydraulic resistence from input vector of coeficients.
    void calc_material_resistance( int type, double coef[6], Iter mat);

};

//#define MT_STOR  4


/// Iterates through materials.
#define FOR_MATERIALS_IT(base,i) \
    for( MaterialDatabase::Iter i = (base).begin(); \
        i != (base).end(); \
        i++)


#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
