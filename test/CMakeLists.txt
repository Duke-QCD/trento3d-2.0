# Create a meta-target which will depend on all the actual testing targets.
# Target 'test' is reserved by CMake.
add_custom_target(trento3dtests)

# Run the compiled executable normally and let the sanitizers do their thing.
# Even if compiled without sanitizers, this won't hurt.
add_custom_target(trento3dsanitizer COMMAND ${PROJECT_NAME} p Pb 10 > /dev/null VERBATIM)
add_dependencies(trento3dtests trento3dsanitizer)

# Download the Catch test header automatically.
set(CATCH_FILE "${CMAKE_CURRENT_SOURCE_DIR}/catch.hpp")
set(SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/download_catch.cmake")
file(WRITE ${SCRIPT_FILE}
"message(STATUS \"Downloading Catch test header\")

file(DOWNLOAD
  \"https://github.com/catchorg/Catch2/releases/download/v2.4.0/catch.hpp\"
  \"${CATCH_FILE}\"
  TIMEOUT 20
  STATUS status
  TLS_VERIFY ON)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR
    \"download failed: code \${status_code} \${status_string}\")
endif()
")
add_custom_command(OUTPUT ${CATCH_FILE} COMMAND ${CMAKE_COMMAND} -P ${SCRIPT_FILE})

# Compile the Catch test executable.
set(TEST_EXE "test-${PROJECT_NAME}")
add_executable(${TEST_EXE} EXCLUDE_FROM_ALL
  catch.hpp
  catch.cxx
  util.cxx
  test_collider.cxx
  test_event.cxx
  test_fast_exp.cxx
  test_nucleon.cxx
  test_nucleus.cxx
  test_output.cxx
)
target_link_libraries(${TEST_EXE} ${LIBRARY_NAME} ${Boost_LIBRARIES} ${HDF5_LIBRARIES})

# Add a target to actually run the tests.
add_custom_target(trento3dcatch COMMAND ${TEST_EXE})
add_dependencies(trento3dtests trento3dcatch)

# Coverage.
find_program(GCOVR gcovr)
if(GCOVR)
  add_custom_target(cov
    COMMAND ${GCOVR}
    --root=${CMAKE_SOURCE_DIR}
    --exclude=${CMAKE_CURRENT_SOURCE_DIR}
    --exclude=${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME}.cxx
    VERBATIM
    DEPENDS trento3dcatch)
  add_dependencies(trento3dtests cov)
endif()
