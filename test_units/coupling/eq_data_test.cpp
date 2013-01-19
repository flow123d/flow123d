/*
 * eq_data_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */

/*
 * TODO:
 * - declaration of BCTypes in EqData classes, - need selection for varion BC types even if the type is given by Discrete field
 *
 *   Potrebuju vytvorit cely strom Input::Type pro FieldBase vracejici ruzne emumy, teoreticky bych se obesel s tim aby to vracelo int
 *   a nemuselo tam byt vice instanci sablony <1,1, XX>, ale jen <1,1,unsigned int>. Pro vraceni hodnot by to stacilo, ale jen bych potreboval aby ruzne instance
 *   te tridy (pro flow::bc_type, transport::bc_typ, reaction:soprption_type) mely i ruzne Input::Type stromy, protoze na jejich konci je pokazde jina Selection
 *   (jinak jsou ty stromy uplne stejne) to by vyzadovalo, nejakou obezlicku v systemu Input::Type (kopirovani celeho podstromu ...)
 *
 *   1) Potomek od Field -> DiscreteField
 *      pretizi metodu get_input_type( Selection ), a vola metody get_input_type pro FieldBase, ta vola podobne metody pro vsechny Field
 *      vsechno dohromady to vytvori samostatnou kopii Input::Type  stromu pro konkretni Selection, takze kazdy klic typu DiscreteField bude mit unikatni strom
 *   2) Pro cteni se pouzije tento unikatni strom, cely abytek kodu zustane v jedne instanci vracejici "int"
 *
 *  TODO:
 *    - stale je problem  patrne s tim, ze musim samostatny strom pro instance Enum valued Fields delat cely znovu a nemonu se v nem odkazovat
 *      na staticke promenne (napr. FieldFormula, FieldPython), behem rc.finish() v bulk_input_type to spadne na
 *      tom, ze najaky pointer p_type nejakeho klice je "uninitialized", to normalne nemuze nastat, takze je mozne ,ze se tam nejak dostal
 *      nejaky spatny kus pameti
 *
 *    - 1) v Input::Types oddelit declare_key s primym a lazy vytvorenim klice, urcite v Record, mozna nejak i v Array (tam problem s tim, ze konstruktor se jmenuje vzdy stejne
 *         a rozliseni jen hlavickou je nedostatecne
 *
 *      2) Field ma pole bool , ktere indikuje, zda jde o Filed vracejici Enum.
 *         Pokud Filed neni Enum, EqDataBase vola  get_input_type - ta vraci
 *         referenci na statickou promennou a v EqDataBase se pak pouzije Record::declare_key( ma parametr predany hodnotou) - prime vytvoreni
 *
 *         Pokud Filed je Enum, tak EqDataBase vola make_input_tree(), a pouziva declare_key_lazy( predava se mu pointer ) pro vytvoreni klice,
 *
 *      3) do Lazy Types davat shared pointery primo na kopie registrovanych objektu, pokusit se to udelat jako soucast TypeBase, predtim
 *         si dat dohromady problemy ktere tim chceme vyresit, problemy: 1) bezpecnost, do Lazy types by se tak dostaly i kopie lokalnich recordu,
 *         cimz by se zabranilo tomu aby se zneplatnily pripadne pointery na ne v recordech a pod. 2) Zjednoduseni - neni treba vytvaret umele recordy jen
 *         s RecordData a pod. 3) LazyTypes implementace ma jedinou rozsahlejsi metodu a to je finish(), jednodussi bude mit v TypeBase pole
 *         pro registarici "lazy_types" (jako shared pointery na prislusne tridy, nicmene konvertovane na shared pointery na TypeBase)
 *
 *
 *
 * - EqData poskytuje metody (generic_input_type, boundary_input_type, bulk_input_type) ktere pro instanci (docasnou)
 *   tridy EqData vrati deklaraci prislusneho stromu typu s tim jsou problemy:
 *   1) prakticky problem, je ze vrsek stromu tvori jen docasne promenne, takze nebude fungovat declare_key
 *   2) pokud bych od jedne equation vytvoril vice instanci, tak budu mit i vice instanci typu pro jejich
 *      eq data klice
 * - Meli bychom mit vsechny rovnice pojmenovane? Jake schema mit pro pojmenovani recordu EqData, ma mit kazda instance konkretni rovnice vlastni
 *   Record (nesmysl)? Tedy v nazvu recordu pouzit nazev tridy
 * - registrace Mesh do potomku FieldBase, myslenka je ze instance potomka FiledBase vzdy slouzi pro vyhodnoceni nad nejakou siti /supersiti
 *   ktera je spojena s rovnici, jednak potrebuju nejaky identifikator site pro kontrolu predavanych ElementAccessoru a pak to potrebuju pro
 *   FieldElementwise.
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
        init_pressure=2.2
      }          
  ],
  bc_data=[
      { rid=10,
        bc_type={TYPE="FieldConstant", value = "dirichlet"},
        bc_pressure={
            TYPE="FieldConstant",
            value=1.23
        }
      } 
  ] 
}
)JSON";









class SomeEquation : public testing::Test, EquationBase {
public:

    class EqData : public EqDataBase {
    public:

        enum BC {
            dirichlet,
            neumann,
            robin,
            total_flux
        };

        static IT::Selection bc_type_selection;


        EqData() : EqDataBase("") {
            ADD_FIELD(init_pressure, "Initial condition as pressure");
            ADD_FIELD(cond_anisothropy, "Anisothropic conductivity tensor.", IT::Default("1.0"));
            ADD_FIELD(bc_type,"Boundary condition type, possible values:");
                      bc_type.set_selection(&bc_type_selection);
            ADD_FIELD(bc_pressure,"Dirichlet BC condition value for pressure.");
            ADD_FIELD(init_conc, "Initial condition for the concentration (vector of size equal to n. components");
                      init_conc.set_n_comp(4);
            // ...
        }
        IT::Array boundary_input_type() {
            static IT::Record rec =  generic_input_type("SomeEquation", "", true)
                                .declare_key("piezo_head", bc_pressure.get_input_type(), "" );
            //rec.finish();
            return IT::Array( rec , 1);
        }
        IT::Array bulk_input_type() {
            static IT::Record rec =  generic_input_type("SomeEquation", "", false)
                                .declare_key("init_piezo_head", init_pressure.get_input_type(), "" );
            //rec.finish();
            return IT::Array( rec , 0);
        }

        Region init_from_input_one_region(Input::Record rec, bool bc_regions) {
            Region region=EqDataBase::init_from_input_one_region(rec, bc_regions);
            if (bc_regions) {
                Input::Iterator<Input::AbstractRecord> field_it = rec.find<Input::AbstractRecord>("piezo_head");
                if (field_it) {
                    bc_pressure(region)->init_from_input(*field_it);
                    //bc_pressure(region)=FieldAddGradient<3, FieldValue<3>::Scalar >(bc_pressure(region), gravity);
                }
            }
        }


        Field<3, FieldValue<3>::Scalar > init_pressure;
        Field<3, FieldValue<3>::TensorFixed > cond_anisothropy;
        Field<3, FieldValue<3>::Scalar > conductivity;
        Field<3, FieldValue<3>::Scalar > water_source_density;
        Field<3, FieldValue<3>::Scalar > storativity;
        Field<3, FieldValue<3>::Vector > init_conc;
        Field<3, FieldValue<3>::Enum > sorption_type;

        BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
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

IT::Selection SomeEquation::EqData::bc_type_selection =
              IT::Selection("EqData_bc_Type")
               .add_value(dirichlet, "dirichlet")
               .add_value(neumann, "neumann")
               .add_value(robin, "robin")
               .add_value(total_flux, "total_flux");


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
    Element el_10(3);
    el_10.region_=Region::db().find_id(10);

    // bulk fields
    {
    ElementAccessor<3> elm(&el_0);
    EXPECT_DOUBLE_EQ(1.1, data.init_pressure.value(p, elm) );

    FieldValue<3>::TensorFixed::return_type value = data.cond_anisothropy.value(p, elm);
    EXPECT_DOUBLE_EQ( 1.0, value.at(0,0) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(0,1) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(0,2) );

    EXPECT_DOUBLE_EQ( 0.0, value.at(1,0) );
    EXPECT_DOUBLE_EQ( 1.0, value.at(1,1) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(1,2) );

    EXPECT_DOUBLE_EQ( 0.0, value.at(2,0) );
    EXPECT_DOUBLE_EQ( 0.0, value.at(2,1) );
    EXPECT_DOUBLE_EQ( 1.0, value.at(2,2) );
    }
    {
    ElementAccessor<3> elm(&el_1);
    EXPECT_DOUBLE_EQ(2.2, data.init_pressure.value(p, elm) );
    }

    //boundary fields
    {
        ElementAccessor<3> elm(&el_10);
        EXPECT_EQ( EqData::dirichlet, data.bc_type.value(p, elm) );
        EXPECT_DOUBLE_EQ(1.23, data.bc_pressure.value(p, elm) );
    }
}
