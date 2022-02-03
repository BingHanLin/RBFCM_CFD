project(test)

include(CTest)

# catch2 test
aux_source_directory(. DIR_SRCS)

add_executable(${PROJECT_NAME} ${DIR_SRCS})

target_link_libraries(${PROJECT_NAME}
    RBFCMCFD
)

install(TARGETS ${PROJECT_NAME}
    DESTINATION . )

add_test(
        NAME
        ${PROJECT_NAME}
        WORKING_DIRECTORY
        ${PROJECT_BINARY_DIR}
        COMMAND
        $<TARGET_FILE:${PROJECT_NAME}>
)



