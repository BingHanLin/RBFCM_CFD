project (THIRD_PARTY)

# Add a library with the above sources
aux_source_directory(${PROJECT_SOURCE_DIR}/pugixml DIR_PUGIXML_SRCS)
add_library(${PROJECT_NAME} 
 ${DIR_PUGIXML_SRCS} )


set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/Eigen
                                 ${PROJECT_SOURCE_DIR}/json
                                 ${PROJECT_SOURCE_DIR}/gmm
                                 ${PROJECT_SOURCE_DIR}/nanoflann
                                 ${PROJECT_SOURCE_DIR}/pugixml
                                 ${PROJECT_SOURCE_DIR}/catch2
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

target_include_directories( ${PROJECT_NAME}
INTERFACE ${PROJECT_SOURCE_DIR}/Eigen
           ${PROJECT_SOURCE_DIR}/json
           ${PROJECT_SOURCE_DIR}/gmm
           ${PROJECT_SOURCE_DIR}/nanoflann
           ${PROJECT_SOURCE_DIR}/pugixml
           ${PROJECT_SOURCE_DIR}/catch2
)
