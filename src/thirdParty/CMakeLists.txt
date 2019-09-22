project (THIRD_PARTY)
add_library(${PROJECT_NAME} INTERFACE)

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/Eigen
                                 ${PROJECT_SOURCE_DIR}/json
                                 ${PROJECT_SOURCE_DIR}/gmm
                                 ${PROJECT_SOURCE_DIR}/nanoflann
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

target_include_directories( ${PROJECT_NAME}
INTERFACE ${PROJECT_SOURCE_DIR}/Eigen
           ${PROJECT_SOURCE_DIR}/json
           ${PROJECT_SOURCE_DIR}/gmm
           ${PROJECT_SOURCE_DIR}/nanoflann
)
