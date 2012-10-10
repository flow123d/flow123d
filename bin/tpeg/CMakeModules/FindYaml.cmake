find_path(YAML_INCLUDE_DIR yaml.h /usr/include /usr/local/include )

find_library(YAML_LIBRARY NAMES yaml PATH /usr/lib /usr/local/lib)

if (YAML_INCLUDE_DIR AND YAML_LIBRARY)
  set(YAML_FOUND TRUE)
endif (YAML_INCLUDE_DIR AND YAML_LIBRARY)

if (YAML_FOUND)
  if (NOT Yaml_FIND_QUIETLY)
    message(STATUS "Found Yaml: ${YAML_LIBRARY}")
  endif (NOT Yaml_FIND_QUIETLY)
else (YAML_FOUND)
  if (Yaml_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find Yaml")
  endif (Yaml_FIND_REQUIRED)
endif (YAML_FOUND)
