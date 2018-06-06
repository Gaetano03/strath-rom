add_executable (prove prove.cpp)
target_link_libraries (prove strath-rom ${MANDATORY_LIBRARIES})

add_executable (Main_Vortex_detection Main_Vortex_detection.cpp)
target_link_libraries (Main_Vortex_detection strath-rom ${MANDATORY_LIBRARIES})

add_executable (RBM_flow_3D RBM_flow_3D.cpp)
target_link_libraries (RBM_flow_3D strath-rom ${MANDATORY_LIBRARIES})

add_executable (RBM_flow_update RBM_flow_update.cpp)
target_link_libraries (RBM_flow_update strath-rom ${MANDATORY_LIBRARIES})