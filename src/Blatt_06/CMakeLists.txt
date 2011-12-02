add_executable(exercise_06_3 exercise_06_3.cpp)
target_link_libraries(exercise_06_3 ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

#running executable and gnuplot
#ADD_CUSTOM_COMMAND( TARGET exercise_01_4
#                    POST_BUILD
#                    COMMAND "exercise_01_4"
#                    COMMAND ${GNUPLOT_EXECUTABLE} "plotscript.p"
#                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
#                    COMMENT "Running executable and printing results with gnuplot...")
                    
if(PDFLATEX_COMPILER)
    SET(Blatt_06_Images
#        exercise_01_4.ps
    )
    
    ADD_LATEX_DOCUMENT(Blatt_06.tex IMAGES ${Blatt_06_Images} DEFAULT_PDF MANGLE_TARGET_NAMES)
endif(PDFLATEX_COMPILER)