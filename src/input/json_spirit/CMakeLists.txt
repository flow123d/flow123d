SET(JSON_SPIRIT_SRCS
json_spirit_reader.cpp 
json_spirit_value.cpp 
json_spirit_writer.cpp
)
 
ADD_LIBRARY(json_spirit STATIC ${JSON_SPIRIT_SRCS})

