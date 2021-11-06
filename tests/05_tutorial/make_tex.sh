#!/bin/bash


# first add blank line between MarkDown description and YAML code
# then add indentation to YAML lines (to be produced as verbatim)
# and remove leading ## from MarkDown lines
sed '/^##.*$/{
      $!{ N        # append the next line when not on the last line
        /^##.*\n##/!s/^\(##.*\)\n/\1\n\n/
                   # now test for a successful substitution, otherwise
                   #+  unpaired "a test" lines would be mis-handled
        t sub-yes  # branch_on_substitute (goto label :sub-yes)
        :sub-not   # a label (not essential; here to self document)
                   # if no substituion, print only the first line
        P          # pattern_first_line_print
        D          # pattern_ltrunc(line+nl)_top/cycle
        :sub-yes   # a label (the goto target of the 't' branch)
                   # fall through to final auto-pattern_print (2 lines)
       }
     }' $1 | sed '/^##/!d;/^##/s/^##\(.*\)/\1/' | sed 's/^#/##/' | sed "1d;2s/^title: .*\`\`\(.*\)''\$/# \1\\n\\nFile: \`\`\`$1\`\`\`\\n/;3d" | pandoc -o ${1%.yaml}.tex --toc -N -V geometry="margin=1in" -f markdown-auto_identifiers --filter pandoc-fignos --filter pandoc-tablenos
