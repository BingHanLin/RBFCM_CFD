project(test)

# catch2 test
aux_source_directory(. DIR_SRCS)

add_executable(${PROJECT_NAME} ${DIR_SRCS})

target_link_libraries(${PROJECT_NAME}
    RBFCMCFD
)

install(TARGETS ${PROJECT_NAME}
    DESTINATION . )




