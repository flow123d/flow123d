#!/bin/bash
FILE=${1:-*.c *.cpp *.h *.hpp *.cc *.hh}
ASO=
# K/R style
ASO="$ASO --style=kr"
# convert tabs to spaces
ASO="$ASO --convert-tabs"	
# do not break one line blocks
ASO="$ASO -- one-line=keep-blocks"
# Insert empty lines around unrelated blocks, labels, classes, ...
#ASO=$ASO --break-blocks

ARTISTIC_STYLE_OPTIONS=$ASO

for f in $FILE;
do 
	if test -f $f 
	then
		astyle $f
	fi
done