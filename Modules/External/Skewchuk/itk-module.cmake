
# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component( MY_CURENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file( READ "${MY_CURENT_DIR}/README" DOCUMENTATION )
 
# define the dependencies of the include module and the tests
itk_module(ITKSKEWCHUK
  DESCRIPTION
    "${DOCUMENTATION}"
)
