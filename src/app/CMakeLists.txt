project(RBFCMCFDAPP)

aux_source_directory(. DIR_MAIN_SRCS)

add_executable(${PROJECT_NAME} ${DIR_MAIN_SRCS})

target_link_libraries(${PROJECT_NAME}
    RBFCMCFD
)

install(TARGETS ${PROJECT_NAME}
    DESTINATION . )


