#  these are defined in findLATEX.cmake:
#
#  LATEX_COMPILER:       path to the LaTeX compiler
#  PDFLATEX_COMPILER:    path to the PdfLaTeX compiler
#  BIBTEX_COMPILER:      path to the BibTeX compiler
#  MAKEINDEX_COMPILER:   path to the MakeIndex compiler
#  DVIPS_CONVERTER:      path to the DVIPS converter
#  PS2PDF_CONVERTER:     path to the PS2PDF converter
#  LATEX2HTML_CONVERTER: path to the LaTeX2Html converter

find_package(LATEX)

if    (NOT LATEX_COMPILER)
  MESSAGE(STATUS "Latex command LATEX_COMPILER not found.")
  set(no_latex)
endif (NOT LATEX_COMPILER)
if    (NOT PDFLATEX_COMPILER)
  MESSAGE(STATUS "Latex command PDFLATEX_COMPILER not found.")
  set(no_latex)
endif (NOT PDFLATEX_COMPILER)
if    (NOT BIBTEX_COMPILER)
  MESSAGE(STATUS "Latex command BIBTEX_COMPILER not found.")
  set(no_latex)
endif (NOT BIBTEX_COMPILER)

if(no_latex)
  MESSAGE(STATUS "Coud not find LaTeX or some of its commands.")
else()
  # if all the LaTeX commands we need where found we can make the reference manual
  
  set(doc_name "flow123d_doc")

  # Here add every new .tex (or other .txt, .sty etc.) file that you create and include
  # with \input \include \usepackage commands in the documentation.
  set(doc_inputs 
       flow_doc.sty
       convection.tex
       darcy_flow.tex
       decay.tex
  #    flow_ini.tex
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
endif(no_latex)
