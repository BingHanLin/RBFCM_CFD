project (RBFCM_CFD)

# Add a library with the above sources
aux_source_directory(${PROJECT_SOURCE_DIR}/geometry/src DIR_MAIN_SRCS)

add_library(${PROJECT_NAME} ${DIR_MAIN_SRCS})

include_directories(${THIRD_PARTY_INCLUDE_DIRS})

target_include_directories( ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/geometry/include
)
