
# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component( MY_CURENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file( READ "${MY_CURENT_DIR}/README" DOCUMENTATION )

# itk_module() defines the module dependencies in ITKPointInCircle
# ITKPointInCircle depends on ITKCommon
# The testing module in ITKPointInCircle depends on ITKTestKernel
# and ITKMetaIO(besides ITKPointInCircle and ITKCore)
 
# define the dependencies of the include module and the tests
itk_module(ITKSkewchuk
DEPENDS
    ITKCommon
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
  DESCRIPTION
    "${DOCUMENTATION}"
)
