add_executable(rentiere rentiere.cc)
target_link_libraries(rentiere ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

#running executable and gnuplot
ADD_CUSTOM_COMMAND( TARGET rentiere
                    POST_BUILD
                    COMMAND "rentiere"
                    COMMAND ${GNUPLOT_EXECUTABLE} "rentiere_plotscript.p"
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    COMMENT "Running executable and printing results with gnuplot...")
                    
add_executable(praediktorkorrektor praediktorkorrektor.cc)
target_link_libraries(praediktorkorrektor ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

if(PDFLATEX_COMPILER)
    SET(ReindeerImages reindeerpopulation.pdf)
    ADD_LATEX_DOCUMENT(Weihnachtsblatt.tex IMAGES ${ReindeerImages} DEFAULT_PDF MANGLE_TARGET_NAMES)
endif(PDFLATEX_COMPILER)
