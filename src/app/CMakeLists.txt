project(MYAPP)

aux_source_directory(. DIR_MAIN_SRCS)

add_executable(${PROJECT_NAME} ${DIR_MAIN_SRCS})

target_link_libraries(${PROJECT_NAME}
RBFCM_CFD
THIRD_PARTY
)


if(CATCH2TEST)
    # catch2 test
    aux_source_directory(./test DIR_TEST_SRCS)

    add_executable(${PROJECT_NAME}_TEST ${DIR_TEST_SRCS})

    target_link_libraries(${PROJECT_NAME}_TEST
    RBFCM_CFD
    THIRD_PARTY
    )
endif()




