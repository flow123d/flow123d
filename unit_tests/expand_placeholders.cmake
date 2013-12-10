#
# In file ${INPUT} expand every line with placeholder $PLACEHOLDER_<name>$ replacing
# placeholder itera
#
# Use template file ${INPUT} to produce ${OUTPUT}. File ${INPUT} contains placeholders in format: $PLACEHOLDER_<filename>$
# One line can contain just one placeholder. Such line is replaced by several lines according to number of lines in file '<filename>'
# where value of the placeholder iterates through lines of that file.
#
# Files with expansion values should be in working directory.
#
#  CONSTRAIN: file ${INPUT} can not contain semi-colons ';'

# Algorithm:
# copy file to temporary
# find all placeholders
# for every placeholder
#   get list from appropriate file
#   for every list item substitute every line

file(READ "${INPUT}" input_str)

# get placeholders
string(REGEX MATCHALL "\\$PLACEHOLDER_[-a-zA-Z0-9_.]*\\$" placeholders "${input_str}")
list(REMOVE_DUPLICATES placeholders)

foreach(placeholder ${placeholders})
  message("PH: ${placeholder}")
  string(REGEX REPLACE "\\$PLACEHOLDER_([-a-zA-Z0-9_.]*)\\$" "\\1" filename "${placeholder}")
  
  string(REGEX MATCH "[^\n]*\\$PLACEHOLDER_${filename}\\$[^\n]*\n" match "${input_str}")
  message("match: ${match}")
  
  message("FN: ${filename}")
  file(READ "${filename}" values)
  STRING(REPLACE "\n" ";" values "${values}")
  foreach(value ${values})
    message("VAL: ${value}")
    # duplicate lines containing placeholder and in the first of two new lines replace placeholder with value
    string(REGEX REPLACE "([^\n]*)(\\$PLACEHOLDER_${filename}\\$)([^\n]*\n)" "\\1${value}\\3\\1\\2\\3" input_str "${input_str}")
  endforeach()
  # remove lines containing placeholder
  string(REGEX REPLACE "([^\n]*)(\\$PLACEHOLDER_${filename}\\$)([^\n]*\n)" "" input_str "${input_str}")
endforeach()

file(WRITE "${OUTPUT}" "${input_str}")
