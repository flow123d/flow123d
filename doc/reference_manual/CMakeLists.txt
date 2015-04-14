# Separate project for making reference manual PDF from 
# Latex sources. UseLatex just sends errors when some
# necessary tool is missing. This wrapper project can 
# be called through execute_process
project(Flow123d_Latex_doc)
cmake_minimum_required(VERSION 2.8)

  set(doc_name "flow123d_doc")

  # Here add every new .tex (or other .txt, .sty etc.) file that you create and include
  # with \input \include \usepackage commands in the documentation.
  set(doc_inputs 
       flow_doc.sty
       ./packages/etoolbox.sty
       convection.tex
       darcy_flow.tex
       chemistry.tex
       dg.tex
  #    flow_ini.tex
       flow_version.tex
       heat.tex
       input_files.tex
       input_reference.tex
       JSON_input.tex
       lumped_mh.tex
       numerical.tex
       output.tex
       pos_12.tex
       reaction_term_numerical.tex
       semchem.tex
       reaction_term.tex
       tests.tex
       transport_model.tex
       tutorial.tex
       version_comparision.tex
       
       # Also add latex table or image if you have it in a separate file.
       tests_graphics/bench_mel_table.tex
  )

  set(doc_img_dirs
      tests_graphics 
      figures
  )

  MESSAGE(STATUS "Targets for LaTeX documentation...")
  INCLUDE(${CMAKE_SOURCE_DIR}/../../CMake/Modules/UseLATEX.cmake)
  MESSAGE(STATUS "UseLATEX.cmake included.")

  ADD_LATEX_DOCUMENT(${doc_name}.tex
		   INPUTS ${doc_inputs}
		   BIBFILES ${doc_name}.bib
		   IMAGE_DIRS ${doc_img_dirs}
		   DEFAULT_PDF
  )
  