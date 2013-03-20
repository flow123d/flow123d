#!/bin/bash
#
#  add_doc_replace.sh <rules> <file_in> [ <file_out> ]
#
# Every line in <rules> has form \AddDoc{...}<text to substitute>. This is a rule to replace
# '\AddDoc{...}' with <text to substitute>,
#
# If both <file_in> and <file_out> are specified,  we apply all rules to the <file_in>
# and write the result into <file_out>.
#
# If only <file_in> is specified, we update the <rules> file according to the list in <file_in> which
# should consist of '\AddDoc{..}' patterns on every line. The <rule> file is applied in different way, 
# every pattern in <file_in> is replaced by the line form <rules> file that contains the same pattern.
# That means <file_in> is the new list of patterns, but rules are reused where it is possible.

#set -x
# file_in
RULES="$1"
RULES_DOUBLE_BS="$1.double"
cat "$RULES" | sed 's/\\/\\\\/g' > "$RULES_DOUBLE_BS"

if [ -z "$3" ] 
then
  # update
  OUT="$RULES"
  UPDATE="yes"
else
  # replace
  OUT="$3"
fi
IN="$2"
OUTTMP="$OUT.tmp"
cp "${IN}" "${OUT}"


while read -r line  
do
  #echo "${line}"
  # parse line
  pattern="${line%%"}"*}}"
  replace="${line#*"}"}"
  if [ -n "$UPDATE" ] 
  then replace="${line}"
  fi  
  #trim spaces
  read -rd '' replace <<< "${replace}"

  # double backslash
  #pattern="${pattern/\\/\\\\}"
  cat "${OUT}" | sed "s/${pattern}/${replace}/" > "${OUTTMP}"  
  mv "${OUTTMP}" "${OUT}"
done  < "${RULES_DOUBLE_BS}"

