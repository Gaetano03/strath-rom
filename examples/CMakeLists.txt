add_executable (prove prove.cpp)
target_link_libraries (prove strath-rom ${MANDATORY_LIBRARIES})

add_executable (Main_Vortex_detection Main_Vortex_detection.cpp)
target_link_libraries (Main_Vortex_detection strath-rom ${MANDATORY_LIBRARIES})

add_executable (RBM_surface_flow RBM_surface_flow.cpp)
target_link_libraries (RBM_surface_flow strath-rom ${MANDATORY_LIBRARIES})

add_executable (RBM_flow RBM_flow.cpp)
target_link_libraries (RBM_flow strath-rom ${MANDATORY_LIBRARIES})