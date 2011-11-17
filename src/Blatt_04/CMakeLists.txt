add_executable(exercise_03_3 exercise_03_3.cpp)
target_link_libraries(exercise_03_3 ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

#running executable and gnuplot
ADD_CUSTOM_COMMAND( TARGET exercise_03_3
                    POST_BUILD
                    COMMAND "exercise_03_3"
                    COMMAND ${GNUPLOT_EXECUTABLE} "plotscript.p"
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    COMMENT "Running executable and printing results with gnuplot...")