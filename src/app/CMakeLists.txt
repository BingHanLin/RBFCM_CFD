project(MYAPP)

aux_source_directory(. DIR_MAIN_SRCS)
add_executable(${PROJECT_NAME} ${DIR_MAIN_SRCS})

target_link_libraries(${PROJECT_NAME}
RBFCM_CFD
THIRD_PARTY
)

