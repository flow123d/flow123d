# Separate project for making reference manual PDF from 
# Latex sources. UseLatex just sends errors when some
# necessary tool is missing. This wrapper project can 
# be called through execute_process
project(Flow123d_Latex_doc)
cmake_minimum_required(VERSION 2.8)

  set(doc_name "flow123d_doc")

  # assume the TUTORIALS variable with list of tutorial tex files.
  # call e.g.: cmake -D TUTORIALS="1.tex 2.tex" . 

  # Here add every new .tex (or other .txt, .sty etc.) file that you create and include
  # with \input \include \usepackage commands in the documentation.
  string(REPLACE " " ";" TUTORIALS_LIST ${TUTORIALS})
  set(doc_inputs 
       flow_doc.sty
       ./packages/etoolbox.sty
       ./packages/multirow.sty
       ${TUTORIALS_LIST}
       abstract_models.tex
       convection.tex
       darcy_flow.tex
       chemistry.tex
       dg.tex
       flow_version.tex
       getting_started.tex
       heat.tex
       input_files.tex
       input_reference.tex
       JSON_input.tex
       lumped_mh.tex
       mechanics.tex
       numerical.tex
       output.tex
       reaction_term_numerical.tex
       semchem.tex
       reaction_term.tex
       tests.tex
       transport_model.tex
       tutorial.tex
       tutorial_new.tex
       version_comparision.tex
       
       # Also add latex table or image if you have it in a separate file.
       tests_graphics/bench_mel_table.tex
  )
  message(STATUS ${doc_inputs})
  set(doc_img_dirs
      tests_graphics 
      figures
      tutor_figures
  )
  
  
#   # generate latex files from tutorials in yaml format
#   set(tutorial_files
#        01_column.yaml
#        02_column_transport.yaml
#        03_tunnel.yaml
#        04_frac_diffusion.yaml
#        05_frac_sorption.yaml
#        06_frac_dualpor.yaml
#        07_heat.yaml
#   )
#   add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/make_tex.sh
#      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/05_tutorial/make_tex.sh
#      COMMAND "cp" make_tex.sh ${CMAKE_CURRENT_BINARY_DIR}
#      COMMAND "cp" -r images ${CMAKE_CURRENT_BINARY_DIR}
#      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/05_tutorial/
#   )
#   foreach(f ${tutorial_files})
# #    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f}
# #       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/05_tutorial/${f}
# #       COMMAND "cp" ${f} ${CMAKE_CURRENT_BINARY_DIR}
# #       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/05_tutorial/
# #    )
# #    add_custom_command(OUTPUT ${f}.tex
# #       DEPENDS ${f} make_tex.sh
# #       COMMAND "${CMAKE_CURRENT_BINARY_DIR}/make_tex.sh" ${f}
# #       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
# #    )
#     add_custom_command(OUTPUT ${f}.tex
#        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/05_tutorial/${f} make_tex.sh
#        COMMAND "cp" ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/05_tutorial/${f} ${CMAKE_CURRENT_BINARY_DIR}
#        COMMAND "${CMAKE_CURRENT_BINARY_DIR}/make_tex.sh" ${f}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     )
# 
#     list(APPEND doc_inputs ${f}.tex)
#   endforeach(f)

  
  

  MESSAGE(STATUS "Targets for LaTeX documentation...")
  INCLUDE(${CMAKE_SOURCE_DIR}/../../CMake/Modules/UseLATEX.cmake)
  MESSAGE(STATUS "UseLATEX.cmake included.")

  ADD_LATEX_DOCUMENT(${doc_name}.tex
		   INPUTS ${doc_inputs}
		   BIBFILES ${doc_name}.bib
		   IMAGE_DIRS ${doc_img_dirs}
		   DEFAULT_PDF
  )
  
