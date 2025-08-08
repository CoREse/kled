set(TEST_RESOURCES
    call_test.sh
    call_test.fa.gz
    call_test.cram
    call_test.cram.crai
    call_test.answer.vcf
)

foreach(resource IN LISTS TEST_RESOURCES)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${resource} AND NOT IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${resource})
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/${resource}
            ${CMAKE_CURRENT_BINARY_DIR}/${resource}
            COPYONLY
        )
    elseif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${resource})
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${resource}/
             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${resource}
        )
    endif()
endforeach()

add_test(
  NAME CallTest
  COMMAND bash call_test.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)