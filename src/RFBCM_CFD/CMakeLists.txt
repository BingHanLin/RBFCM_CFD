project (RBFCM_CFD)

# Add a library with the above sources
aux_source_directory(${PROJECT_SOURCE_DIR}/mesh/src DIR_MESH_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/boundaryConditions/src DIR_BOUNDARYCONDITIONS_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/rbfBasis/src DIR_RBFBASIS_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/simulationDomain/src DIR_SIMULATIONDOMAIN_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/controlData/src DIR_CONTROLDATA_SRCS)

add_library( ${PROJECT_NAME}
    ${DIR_MESH_SRCS} 
    ${DIR_BOUNDARYCONDITIONS_SRCS} 
    ${DIR_RBFBASIS_SRCS} 
    ${DIR_SIMULATIONDOMAIN_SRCS}
    ${DIR_CONTROLDATA_SRCS}
)

include_directories(${THIRD_PARTY_INCLUDE_DIRS})

target_include_directories( ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/mesh/include
    ${PROJECT_SOURCE_DIR}/boundaryConditions/include
    ${PROJECT_SOURCE_DIR}/rbfBasis/include
    ${PROJECT_SOURCE_DIR}/simulationDomain/include
    ${PROJECT_SOURCE_DIR}/controlData/include
    ${PROJECT_SOURCE_DIR}/enums/include
    ${PROJECT_SOURCE_DIR}/utility/include

)
