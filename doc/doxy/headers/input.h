/**
 * @defgroup input_mod  Input Classes
 *
 * \section classes Overview of data input classes
 * We assume that the input data of the program can be represented by a tree.
 * The leave nodes of the tree contain actual data and has one of the 'scalar' types
 * represented by the classes in Input::Type namespace : Bool, Integer, Double, String, 
 * and FileName. The branching nodes of the input tree 
 * can be either of type Array or type Record (or Abstract).
 * Classes derived from \p Input::Type::TypeBase only describes structure of the input tree
 * which is independent of actual format of the input data. 
 * All nodes of tree allows definition of input attributes, which are optional for every node.
 * Base attributes are defined in Input::Type::Attributes, special attributes of Flow123D in
 * FlowAttributes. These structures contain common attributes, user can define other else.
 * 
 * Instances of the classes from 
 * the Input::Type namespace form a \e declaration-tree  that is later used
 * by a reader of the particular file format (currently we support YAML and JSON through
 * the Input::ReaderToStorage reader) to interpret the input data, check its structure, possibly use default values
 * and put the raw data into an intermediate \e storage-tree formed be Input:Storage classes. 
 * 
 * Finally, 
 * the data are accessible through accessors Input::Record, Input::Array, and Input::AbstractRecord.
 * These accessors holds pointers into declaration tree as well as into the data storage tree 
 * and provides unified access to the data. 
 * 
 * Furthermore, the \e declaration-tree can output itself provided a basic documentation of the 
 * input data structure, that is consistent with the actual state.
 * 
 * Here is simple scheme of information exchange between classes:
 * 
 * @dot
 digraph G
{
  graph[rankdir="LR",bgcolor="transparent"];
  
  edge [fontname="FreeSans",fontsize=15,labelfontname="FreeSans",labelfontsize=10];
  node [fontname="FreeSans",fontsize=15,
        shape=record,height=0.2,width=0.4,
        color="black", fillcolor="white", style="filled"];

  DeclarationTree [label="Declaration tree\n(Input::Type...)",URL="\ref input_types"];
  InputFile [label="Input file\n(YAML, JSON)"];
  Reader [label="Particular reader\n(Input::ReaderToStorage)",URL="\ref Input::ReaderToStorage"];
  Storage [label="Data storage tree\n(Input::StorageBase...)",URL="\ref Input::StorageBase"];
  Accessors [label="Data accessors\n(Input::Record, Input::Array, ...)",URL="\ref input_accessors"];
  Doc [label="Input documentation"];
  DataUsage [label="Data usage"];

  {rank=same; InputFile Reader Storage Accessors DataUsage}
  {rank=same; DeclarationTree Doc }
  
  DeclarationTree -> Reader [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  DeclarationTree -> Accessors [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  DeclarationTree -> Doc [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  InputFile -> Reader [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  Reader -> Storage [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  Storage -> Accessors [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  Accessors -> DataUsage [color="black",fontsize=10,style="solid",fontname="FreeSans"];
}
 * @enddot
 * 
 * \section practical_usage Practical usage
 * 
 * In order to use input interface you have to create \e declaration-tree, namely Selection, Record, and Abstract types
 * has to be declared by specification of its keys. We suggest to have a static method that returns Input::Type::Record 
 * for ever class which is initialized from input interface. Example of usage follows:
 * 
 @code
    class Mesh {
      static Input::Type::Record get_input_type() {
          using namespace Input::Type; 
          return Record("Mesh", "Full description of computational mesh.")
              .declare_key("mesh_name", String(), Default::optional(),"Optional name of the mesh")
              .declare_key("mesh_file", FileName::input(), Default("mesh.msh"), "Principial mesh input file")
              .declare_key("materials_to_remove", Array(Integer()), "Removes elements with material ID that appears in this list.")
              .close();
      }  
      
      Mesh( Input::Record input) 
      : name( input.val<string>("mesh_name") )
      {
        string fname = input.val<FilePath>("mesh_file");
        boundary = new Boundary( input.val<Input::Record>("boundary") );
        std::vector<int> ids_to_remove;
        input.val<Input::Array>("materials_to_remove").copy_to( ids_to_remove );
      }  
    }  
 @endcode  
 * 
 * The accessor for the root Input::Record is provided by the reader class ( ReaderToStorage.get_root_interface<>() )
 *
 */
