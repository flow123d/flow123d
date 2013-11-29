
set(doc_name "flow123d_doc")

# Here add every new .tex (or other .txt, .sty etc.) file that you create and include
# with \input \include \usepackage commands in the documentation.
set(doc_inputs 
       flow_doc.sty
       convection.tex
       darcy_flow.tex
       decay.tex
#       flow_ini.tex
       input_files.tex
       input_reference.tex
       JSON_input.tex
       numerical.tex
       output.tex
       pos_12.tex
       semchem.tex
       sorpce.tex
       tests.tex
       transport_model.tex
       version_comparision.tex
       
       # Also add latex table or image if you have it in a separate file.
       tests_graphics/bench_mel_table.tex
)

set(doc_img_dirs
      tests_graphics 
      figures
)

MESSAGE(STATUS "Buiding LaTeX documentation.")
INCLUDE(${PROJECT_SOURCE_DIR}/CMake/Modules/UseLATEX.cmake)
MESSAGE(STATUS "UseLATEX.cmake included.")

ADD_LATEX_DOCUMENT(${doc_name}.tex
		   INPUTS ${doc_inputs}
		   BIBFILES ${doc_name}.bib
		   IMAGE_DIRS ${doc_img_dirs}
		   DEFAULT_PDF
)
# NO_DEFAULT - does not build anything with make all, "make pdf" must be called
# Put images in INPUTS, then no conversion is done
# or convert them to pdf and include the image directory

add_custom_target(ref-doc
  COMMAND make pdf
)

