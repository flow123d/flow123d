
set(doc_name "flow123d_doc")
set(doc_inputs 
       flow_doc.sty
       bcd_10.tex
       convection.tex
       darcy_flow.tex
       decay.tex
       flow_ini.tex
       input_files.tex
       input_reference.tex
       JSON_input.tex
       msh_20.tex
       mtr_10.tex
       ngh_10.tex
       numerical.tex
       output.tex
       pos_12.tex
       semchem.tex
       sorpce.tex
       src_10.tex
       tests.tex
       transport_model.tex
       version_comparision.tex
       
       #table
       tests_graphics/bench_mel_table.tex
 )

MESSAGE(STATUS "Buiding LaTeX documentation.")
INCLUDE(${PROJECT_SOURCE_DIR}/CMake/Modules/UseLATEX.cmake)
MESSAGE(STATUS "UseLATEX.cmake included.")

ADD_LATEX_DOCUMENT(${doc_name}.tex
		   INPUTS ${doc_inputs}
		   BIBFILES ${doc_name}.bib
		   IMAGE_DIRS tests_graphics figures
		   DEFAULT_PDF
    )
    
add_custom_target(doc-pdf
  COMMAND make pdf
)
# NO_DEFAULT - does not build anything, "make pdf" must be called
# put images in INPUTS, then no conversion is done
