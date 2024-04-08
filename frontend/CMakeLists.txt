# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/frontend_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/frontend_config.h"
)

# set binary path and lib storage path for compiled files.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set(HEADER_FILES
        config_values.h
        misa_md.h
        def_config_values.h
        config_parser.h
        device.h
        md_simulation.h
        io/atom_dump.h
        io/atom_dump_types.h
        io/buffered_io.h
        io/output_copy.h
        io/output_base_interface.h
        io/output_dump.h
        )

set(SOURCE_FILES
        config_values.cpp
        config_values_thermodynamic.hpp
        misa_md.cpp
        config_parser.cpp
        main.cpp
        md_simulation.cpp
        io/atom_dump.cpp
        io/atom_dump_types.cpp
        io/buffered_io.cpp
        io/output_copy.cpp
        io/output_dump.cpp
        )

add_library(${MD_FRONTEND_LIB_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(
        ${MD_FRONTEND_LIB_NAME}
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(
        ${MD_FRONTEND_LIB_NAME}
        PUBLIC
        yaml-cpp
        ${MD_LIB_NAME}
)

add_executable(${EXECUTE_BIN_NAME} main.cpp)

target_link_libraries(${EXECUTE_BIN_NAME}
        PRIVATE ${MD_FRONTEND_LIB_NAME} ${EXTRA_LIBS}
        )

# install bin and libs
install(TARGETS ${EXECUTE_BIN_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
