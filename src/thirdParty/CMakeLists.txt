project (thirdParty)

aux_source_directory(${PROJECT_SOURCE_DIR}/pugixml DIR_PUGIXML_SRCS)

add_library(${PROJECT_NAME} STATIC
    ${DIR_PUGIXML_SRCS}
)

target_include_directories(${PROJECT_NAME}
    INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Eigen>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/json>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/gmm>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/nanoflann>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/pugixml>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/catch2>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cxxopts>  
        $<INSTALL_INTERFACE:include>
)