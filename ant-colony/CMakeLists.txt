#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(antColony src/AntColony.cpp)
target_link_libraries(antColony PUBLIC
        ${LIB_GMDS_IG} ${LIB_GMDS_IO} ${LIB_GMDS_IG_ALGO})
target_compile_features(antColony PUBLIC cxx_std_14)

