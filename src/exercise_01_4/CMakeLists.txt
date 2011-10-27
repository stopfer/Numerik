add_executable(exercise_01_4 exercise_01_4.cpp)
target_link_libraries(exercise_01_4 ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

#running executable and gnuplot
ADD_CUSTOM_COMMAND( TARGET exercise_01_4
                    POST_BUILD
                    COMMAND "exercise_01_4"
                    COMMAND ${GNUPLOT_EXECUTABLE} "plotscript.p"
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    COMMENT "Running executable and printing results with gnuplot...")