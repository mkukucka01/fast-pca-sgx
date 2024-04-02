# Copyright (c) Open Enclave SDK contributors.
# Licensed under the MIT License.
add_custom_command(
  OUTPUT fast_pca_u.h fast_pca_u.cpp fast_pca_args.h
  DEPENDS ${CMAKE_SOURCE_DIR}/fast_pca.edl
  COMMAND
    openenclave::oeedger8r --untrusted ${CMAKE_SOURCE_DIR}/fast_pca.edl
    --search-path ${OE_INCLUDEDIR} --search-path
    ${OE_INCLUDEDIR}/openenclave/edl/sgx)

add_executable(fast_pca_host host.cpp
                               ${CMAKE_CURRENT_BINARY_DIR}/fast_pca_u.cpp)

if (WIN32)
  copy_oedebugrt_target(fast_pca_host_oedebugrt)
  add_dependencies(fast_pca_host fast_pca_host_oedebugrt)
endif ()

target_include_directories(
  fast_pca_host PRIVATE # Needed for the generated file fast_pca_u.h
                          ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(fast_pca_host openenclave::oehost)