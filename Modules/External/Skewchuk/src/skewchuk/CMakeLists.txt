project(ITK3P_SKEWCHUK)

include_regular_expression("^(itk|predicates|skewchuk).*$")

include_directories(BEFORE
               "${CMAKE_CURRENT_BINARY_DIR}"
               "${ITK3P_SKEWCHUK_SOURCE_DIR}")

# memmgr back ends: compile only one of these into a working library
# (For now, let's use the mode that requires the image fit into memory.
# This is the recommended mode for Win32 anyway.)
set(systemdependent_SRCS jmemnobs.c)

# library object files common to compression and decompression
set(skewchuk_SRCS
predicates.c
)

#if(WIN32)
#  if(BUILD_SHARED_LIBS)
#    add_definitions(-DSKEWCHUKDLL)
#  else(BUILD_SHARED_LIBS)
#    set(SKEWCHUKSTATIC 1)
#  endif(BUILD_SHARED_LIBS)
#endif(WIN32)

configure_file(${ITK3P_SKEWCHUK_SOURCE_DIR}/skewchukDllConfig.h.in
               ${ITK3P_SKEWCHUK_BINARY_DIR}/skewchukDllConfig.h)

add_library(itkskewchuk ${skewchuk_SRCS})

# Apply user-defined properties to the library target.
if(ITK_LIBRARY_PROPERTIES)
  set_target_properties(itkskewchuk PROPERTIES ${ITK_LIBRARY_PROPERTIES})
endif(ITK_LIBRARY_PROPERTIES)

if(NOT ITK_INSTALL_NO_LIBRARIES)
  install(TARGETS itkskewchuk
    EXPORT ${ITK3P_INSTALL_EXPORT_NAME}
    RUNTIME DESTINATION ${ITK3P_INSTALL_RUNTIME_DIR} COMPONENT RuntimeLibraries
    LIBRARY DESTINATION ${ITK3P_INSTALL_LIBRARY_DIR} COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION ${ITK3P_INSTALL_ARCHIVE_DIR} COMPONENT Development
    )
endif(NOT ITK_INSTALL_NO_LIBRARIES)

if(NOT ITK_INSTALL_NO_DEVELOPMENT)
  install(FILES
    ${ITK3P_SKEWCHUK_BINARY_DIR}/skewchukDllConfig.h
    ${ITK3P_SKEWCHUK_SOURCE_DIR}/predicates.c
    DESTINATION ${ITK3P_INSTALL_INCLUDE_DIR}/itkskewchuk
    COMPONENT Development)
endif(NOT ITK_INSTALL_NO_DEVELOPMENT)
