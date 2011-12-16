add_executable(exercise_7 exercise_7.cpp)
target_link_libraries(exercise_7 ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})


if(PDFLATEX_COMPILER)
    ADD_LATEX_DOCUMENT(Blatt_07.tex  DEFAULT_PDF MANGLE_TARGET_NAMES)
endif(PDFLATEX_COMPILER)
