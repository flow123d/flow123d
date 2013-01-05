/*
 * eq_data_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */

/*
 * TODO:
 * - declaration of BCTypes in EqData classes, - need selection for varion BC types even if the type is given by Discrete field
 * - EqData poskytuje metody (generic_input_type, boundary_input_type, bulk_input_type) ktere pro instanci (docanou)
 *   tridy EqData vrati deklaraci prislusneho stromu typu s tim jsou problemy:
 *   1) prakticky problem, je ze vrsek stromu tvori jen docasne promenne, takze nebude fungovat declare_key
 *   2) pokud bych od jedne equation vytvoril vice instanci, tak budu mit i vice instanci typu pro jejich
 *      eq data klice
 * - Meli bychom mit vsechny rovnice pojmenovane? Jake schema mit pro pojmenovani recordu EqData, ma mit kazda instance konkretni rovnice vlastni
 *   Record (nesmysl)? Tedy v nazvu recordu pouzit nazev tridy
 */

#include <gtest/gtest.h>
#include <vector>
#include <boost/foreach.hpp>

#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"

#include "fields/field_base.hh"
#include "coupling/equation.hh"
#include "mesh/region.hh"

using namespace std;
namespace IT=Input::Type;

// Test data strings
const string eq_data_input = R"JSON(
{ 
  bulk_data=[
      { rid=0,
          init_pressure={
              TYPE="FieldConstant",
              value=1.1
          }
      },
      { rid=1,
          init_pressure={
              TYPE="FieldConstant",
              value=2.2
          }
      }          
  ],
  bc_data=[
      { rid=10,
        bc_pressure={
            TYPE="FieldConstant",
            value=1.0
        }
      } 
  ] 
}
)JSON";



class RegionFieldBase {
public:
    RegionFieldBase(bool bc) : bc_(bc) {}

    void set_name(const string & name)
        { name_ =name; }
    void set_desc(const string & desc)
        { desc_=desc; }
    void set_n_comp( unsigned int n_comp)
        { n_comp_=n_comp; }

    virtual IT::AbstractRecord &get_input_type() =0;

    const std::string &name() const
        { return name_;}
    const std::string &desc() const
        { return desc_;}
    bool is_bc() const
        { return bc_;}
    virtual void init_from_input(Region reg, Input::AbstractRecord rec) =0;
protected:
    std::string name_;
    std::string desc_;
    bool bc_;
    unsigned int n_comp_;

};



template<int spacedim, class Value>
class Field : public RegionFieldBase {
public:
    typedef FieldBase<spacedim, Value> Field_Base;

    IT::AbstractRecord &get_input_type()
        { return Field_Base::input_type; }

    /**
     * Default constructor. This allocates an array that attach actual field instance to every region
     * thus we assumes that at the construction time the RegionDatabase contains at least all regions for which
     * the values of the field will be required.
     *
     */
    Field() : RegionFieldBase(false)  {}
    inline  Field_Base * operator() (Region reg)
        { ASSERT_LESS(reg.idx(), this->region_fields.size());
          return this->region_fields[reg.idx()];
        }

    void init_from_input(Region reg, Input::AbstractRecord rec)
    {
        // initialize table if it is empty, we assume that the RegionDB is closed at this moment
        if (region_fields.size() == 0) region_fields.resize(Region::db().size(), NULL);

        if (region_fields[reg.idx()] !=NULL) {
            xprintf(Warn, "Overwriting existing value on region ID=%d. In initialization of the Field: '%s'.\n", reg.id(), this->name().c_str());
            delete region_fields[reg.idx()];
        }
        region_fields[reg.idx()] = Field_Base::function_factory(rec, 0.0, this->n_comp_);
    }
    /**
     * Returns one value in one given point @p on an element given by ElementAccessor @p elm.
     * It returns reference to he actual value in order to avoid temporaries for vector and tensor values.
     *
     * This method just call the later one @p value(Point, ElementAccessor, Value) and drops the FieldResult.
     */
    inline typename Value::return_type &value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm)
    {
        return region_fields[elm.region().idx()]->value(p,elm);
    }

    /**
     * Pure virtual method. At least this has to be implemented by descendants.
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    inline FieldResult value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm, typename Value::return_type &value)
    {
        return region_fields[elm.region().idx()]->value(p,elm,value);
    }

    /**
     * Returns std::vector of scalar values in several points at once. The base class implements
     * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
     * but this overhead can be negligable for more complex fields as Python of Formula.
     */
    inline void value_list(const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list,
                       std::vector<FieldResult> &result_list)
    {
        region_fields[elm.region().idx()]->value_list(point_list,elm, value_list, result_list);
    }


private:
    std::vector<Field_Base *> region_fields;
};



template<int spacedim, class Value>
class BCField : public Field<spacedim, Value> {
public:
    BCField() { this->bc_=true; }
};



class EqDataBase {
public:
    void add_field( RegionFieldBase *field, const string &name, const string &desc) {
        field->set_name( name );
        field->set_desc( desc );
        field_list.push_back(field);
    }
    /**
     * - @p eq_class_name should be name of the particular equation class, the name of bulk data record has form:
     *   'EqName_BulkData' and record for boundary data has name 'EqName_BoundaryData'. However, these names has
     *   only documentation purpose since these records are not descendants of an AbstractRecord.
     * - we do not finish returned records !!
     */
    IT::Record generic_input_type(const string &eq_class_name, const string &desc, bool bc_regions) {
        string rec_name = eq_class_name + (bc_regions ? "_BoundaryData" : "_BulkData");
        IT::Record rec = IT::Record(rec_name, desc)
                         .declare_key("region", IT::String(), "Region label")
                         .declare_key("rid", IT::Integer(0), "Region ID (alternative to region label)" );

        BOOST_FOREACH(RegionFieldBase * field, field_list)
                    if (bc_regions == field->is_bc()) rec.declare_key(field->name(), field->get_input_type(), field->desc() );


        // intentionally we do not call finish here in order to allow adding keys after the generic ones
        // finish should be called at global level through LazyTypes

        return rec;
    }

    void init_from_input(Input::Array bulk_list, Input::Array bc_list) {
        for(Input::Iterator<Input::Record> it=bulk_list.begin<Input::Record>(); it != bulk_list.end(); ++it)
            init_from_input_one_region(*it, false);
        for(Input::Iterator<Input::Record> it=bc_list.begin<Input::Record>(); it != bc_list.end(); ++it)
            init_from_input_one_region(*it, true);
    }

    Region init_from_input_one_region(Input::Record rec, bool bc_regions) {
        Input::Iterator<string> it = rec.find<string>("region");
        Region reg;

        // get the region
        if (it) {
            // try find region by label
            reg = Region::db().find_label(*it);
            if (! reg.is_valid() ) xprintf(UsrErr, "Unknown region with label: '%s'\n", (*it).c_str());
        } else {
            // try find region by ID
            Input::Iterator<unsigned int> id_it = rec.find<unsigned int>("rid");
            reg = Region::db().find_id(*id_it);
            if (! reg.is_valid() ) xprintf(UsrErr, "Unknown region with id: '%d'\n", *id_it);
        }

        // init all fields on this region
        BOOST_FOREACH(RegionFieldBase * field, field_list) {
            if (bc_regions == field->is_bc()) {
                Input::Iterator<Input::AbstractRecord> field_it = rec.find<Input::AbstractRecord>(field->name());
                if (field_it) {
                    field->init_from_input(reg,*field_it);
                }
            }
        }

        return reg;
    }

protected:
    std::vector<RegionFieldBase *> field_list;
    //static IT::Record *generic_input_type;
};

#define ADD_FIELD(name, desc)     add_field(&name, string(#name), string(desc) )





class SomeEquation : public testing::Test, EquationBase {
public:


    class EqData : public EqDataBase {
    public:
        EqData() {
            ADD_FIELD(init_pressure, "Initial condition as pressure");
            ADD_FIELD(cond_anisothropy, "Anisothropic conductivity tensor.");
            ADD_FIELD(bc_type,"Boundary condition type, possible values:");
            ADD_FIELD(bc_pressure,"Dirichlet BC condition value for pressure.");
            ADD_FIELD(init_conc, "Initial condition for the concentration (vector of size equal to n. components");
            init_conc.set_n_comp(4);
            // ...
        }
        IT::Array boundary_input_type() {
            static IT::Record rec =  generic_input_type("SomeEquation", "", true)
                                .declare_key("piezo_head", bc_pressure.get_input_type(), "" );
            return IT::Array( rec , 1);
        }
        IT::Array bulk_input_type() {
            static IT::Record rec =  generic_input_type("SomeEquation", "", false)
                                .declare_key("init_piezo_head", init_pressure.get_input_type(), "" );
            return IT::Array( rec , 0);
        }

        void init_from_input_one_region(Input::Record rec, bool bc_regions) {
            Region region=EqDataBase::init_from_input_one_region(rec, bc_regions);
            Input::Iterator<Input::AbstractRecord> field_it = rec.find<Input::AbstractRecord>("piezo_head");
            if (field_it) {
                bc_pressure(region)->init_from_input(*field_it);
                //bc_pressure(region)=FieldAddGradient<3, FieldValue<3>::Scalar >(bc_pressure(region), gravity);
            }
        }


        Field<3, FieldValue<3>::Scalar > init_pressure;
        Field<3, FieldValue<3>::TensorFixed > cond_anisothropy;
        Field<3, FieldValue<3>::Scalar > conductivity;
        Field<3, FieldValue<3>::Scalar > water_source_density;
        Field<3, FieldValue<3>::Scalar > storativity;
        Field<3, FieldValue<3>::Vector > init_conc;
        Field<3, FieldValue<3>::Discrete > sorption_type;

        BCField<3, FieldValue<3>::Discrete > bc_type;
        BCField<3, FieldValue<3>::Scalar > bc_pressure; // ?? jak pridat moznost zadat piezo_head, coz by melo initializovat pressure
                                                     // na AddGradient(..)
                                                     // jednak jak deklatovat ten klic, dale jak behem cteni inicializovat pressure
                                                     // tj. potrebuju umet pridat dalsi klice a dalsi inicializace i po generickych funkcich
        BCField<3, FieldValue<3>::Scalar > bc_flux;
        BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;

    };

public:
    SomeEquation() : EquationBase(){

    }

    void get_solution_vector(double*&, unsigned int&) {}
    void get_parallel_solution_vector(_p_Vec*&) {}
protected:
    static Input::Type::Record input_type;
    EqData data;


    virtual void SetUp() {
        // read input string
        std::stringstream ss(eq_data_input);
        Input::JSONToStorage reader;
        reader.read_stream( ss, input_type );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        Region::db().add_region(0,"main_volume",3,false);
        Region::db().add_region(1,"upper_layer",3,false);
        Region::db().add_region(10,"top_surface",3,true);
        data.init_from_input( in_rec.val<Input::Array>("bulk_data"), in_rec.val<Input::Array>("bc_data") );
    }
    virtual void TearDown() {
    };



};

IT::Record SomeEquation::input_type=
        IT::Record("SomeEquation","")
        .declare_key("bc_data", SomeEquation::EqData().boundary_input_type(), IT::Default::obligatory(), ""  )
        .declare_key("bulk_data", SomeEquation::EqData().bulk_input_type(), IT::Default::obligatory(), ""  );



TEST_F(SomeEquation, values) {
    Point<3> p;
    p(0)=1.0; p(1)= 2.0; p(2)=3.0;

    Element el_0(3);
    el_0.region_=Region::db().find_id(0);
    Element el_1(3);
    el_1.region_=Region::db().find_id(1);

    {
    ElementAccessor<3> elm(&el_0);
    EXPECT_DOUBLE_EQ(1.1, data.init_pressure.value(p, elm) );
    }
    {
    ElementAccessor<3> elm(&el_1);
    EXPECT_DOUBLE_EQ(2.2, data.init_pressure.value(p, elm) );
    }
}
