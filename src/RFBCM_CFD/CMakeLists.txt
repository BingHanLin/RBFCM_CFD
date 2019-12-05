project (RBFCM_CFD)

# Add a library with the above sources
aux_source_directory(${PROJECT_SOURCE_DIR}/geometry/src DIR_GEOMETRY_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/rbfBasis/src DIR_RBFBASIS_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/simulationDomain/src DIR_RBFBASIS_SRCS)

add_library(${PROJECT_NAME} ${DIR_GEOMETRY_SRCS} ${DIR_RBFBASIS_SRCS})

include_directories(${THIRD_PARTY_INCLUDE_DIRS})

target_include_directories( ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/geometry/include
    ${PROJECT_SOURCE_DIR}/rbfBasis/include
    ${PROJECT_SOURCE_DIR}/simulationDomain/include
    ${PROJECT_SOURCE_DIR}/enums/include
    ${PROJECT_SOURCE_DIR}/utility/include

)
