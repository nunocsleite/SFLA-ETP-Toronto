# project
#-------------------------------
project(exes)

# Set public header list (add your headers and source files here))
#-------------------------------
set(${PROJECT_NAME}_headers )
set(${PROJECT_NAME}_sources TestMOEA.cpp)

# Create binaries
#-------------------------------
add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_headers} ${${PROJECT_NAME}_sources})

# Link
#-------------------------------
target_link_libraries(${PROJECT_NAME} eo es moeo cma eoutils ga)
target_link_libraries(${PROJECT_NAME} MOEAlib)
