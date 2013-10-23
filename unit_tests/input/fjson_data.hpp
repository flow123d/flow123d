/*
 * fjson_data.hpp
 *
 *  Created on: 27.4.2012
 *      Author: jiri
 */

#ifndef FJSON_DATA_HPP_
#define FJSON_DATA_HPP_

#include <string>

// Test data strings
const string flow_mini_json = R"JSON(
{
    "flow_ini_version"  : "1.0",
    "comment"           : "verified by http://json.parser.online.fr/ to be VALID JSON",

    "global" : {
        "problem_type"  : 1,
        "description"   : "test1",
        "save_step"     : 0.1,
        "density_on"    : false,
        "nothing"       : null
    },

    "input" : {
        "file_type"         : 1,
        "mesh"              : "./input/test1.msh"
    },

    "constants" : {
        "g"       :9.81,
        "rho"     :1000
    },

    "sp" : {
        "drfl"              :1e-009,
        "l_size"            :80
    },

    "z_test_weird_array" : [ [0], { "a" : 1 }, 2, {}, [] ]
}
)JSON";

const string flow_json_comment_parser = R"JSON(
# komentar na zacatku, obsahuje humus { " } " : # \ \\ \ \{ \" \} \: \# 
{
# komentar uvnitr, obsahuje humus { " } " : # \ \\ \ \{ \" \} \: \#


    "text0"           : "text",

#viceradkovy dlouhy komentar \
pokracovani komentare \
jeste dalsi pokracovani komentare \
pokracovani s humusem { " } " : # \ \\ \ \{ \" \} \: \# \
pokracovani bez humusu    
    
    "text1"           : "text", # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
     
    "text2"           : "text" # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    , 
    
    "text3"           : # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    "text",
     
    "text4"          # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
     : "text",

    "text5\""           : "text",
    "text6\"\""           : "text",
    "text7"           : "text\"",
    "text8"           : "text\"\"",
    "text9\""           : "text\"",
    "text10#"           : "text",
    "text11"           : "text#",
    "text12#"           : "text#",

    "record0" : { # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
        "subrecord0"  : 1, # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \# 
        
        "subrecord1"  : 1 # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
        ,
         
        "subrecord2"  : # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \# 
        1,
         
        "subrecord3"  # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
        : 1
         
    }, # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    
    "record1" : { } # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    ,
    
    "record2" : # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \# 
    { },
    
    "record3" # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \# 
    : { },


    "z_array0" : [ [0], { "a" : 1 }, 2, {}, [] ], # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    "z_array1" : [ [0], { "a" : 1 }, 2, {}, [] 
    ],
    
    "z_array2" : [ [0], { "a" : 1 }, 2, {}, [ # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    ] ],
    
    "z_array3" : [ [0], { "a" : 1 }, 2, {}, # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
     [] ],
     
    "z_array4" : [ [0], { "a" : 1 }, 2, { # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    }, [] ],
    
    "z_array5" : [ [0], { "a" : # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
     1 }, 2, {}, [] ],
     
    "z_array6" : [ 
    [0], { "a" : 1 }, 2, {}, [] ], # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
    
    "z_array7" :
     [ [0], { "a" : 1 }, 2, {}, [] ], # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
     
    "z_array8"
     : [ [0], { "a" : 1 }, 2, {}, [] ] # komentar s humusem { " } " : # \ \\ \ \{ \" \} \: \#
}
# komentar na konci, obsahuje humus { " } " : # \ \\ \ \{ \" \} \: \#
)JSON";


const string flow_json_colon_eq = R"JSON(
{
    "key1" = "value",
    "key2" : "value"
}
)JSON";

const string flow_json_quotes = R"JSON(
{
    "key1" = "value1",
    "key2" : "value2",
    "key3" = "value3",
    "key4" : "value4",
     key5  = "value5",
     key6  : "value6"
}
)JSON";

const string flow_json_whitespace_separator = R"JSON(
{
    "key1" = "value1",
    "key2" : "value2",
     key3  = "value3",
     key4  : "value4"
    "key5" = "value5"
    "key6" : "value6"
     key7  = "value7"
     key8  : "value8"
}
)JSON";






#endif /* FJSON_DATA_HPP_ */
