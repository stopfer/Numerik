add_executable(exercise_8 lmm.cc)
target_link_libraries(exercise_8 ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})


if(PDFLATEX_COMPILER)
    ADD_LATEX_DOCUMENT(blatt8.tex  DEFAULT_PDF MANGLE_TARGET_NAMES)
endif(PDFLATEX_COMPILER)
