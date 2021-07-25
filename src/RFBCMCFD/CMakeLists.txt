project (RBFCMCFD)

aux_source_directory(${PROJECT_SOURCE_DIR}/mesh/src DIR_MESH_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/initialConditions/src DIR_INITIALCONDITIONS_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/boundaryConditions/src DIR_BOUNDARYCONDITIONS_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/rbfBasis/src DIR_RBFBASIS_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/simulationDomain/src DIR_SIMULATIONDOMAIN_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/controlData/src DIR_CONTROLDATA_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/conditionPool/src DIR_CONDITIONPOOL_SRCS)
aux_source_directory(${PROJECT_SOURCE_DIR}/utility/src DIR_UTILITY_SRCS)


add_library( ${PROJECT_NAME}
    ${DIR_MESH_SRCS} 
    ${DIR_INITIALCONDITIONS_SRCS}
    ${DIR_BOUNDARYCONDITIONS_SRCS} 
    ${DIR_RBFBASIS_SRCS} 
    ${DIR_SIMULATIONDOMAIN_SRCS}
    ${DIR_CONTROLDATA_SRCS}
    ${DIR_CONDITIONPOOL_SRCS}
    ${DIR_UTILITY_SRCS}
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/mesh/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/initialConditions/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/boundaryConditions/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/rbfBasis/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/simulationDomain/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/controlData/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/enums/include>  
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utility/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/conditionPool/include>  
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME}
    thirdParty
)