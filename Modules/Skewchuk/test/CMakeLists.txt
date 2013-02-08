itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set( ITK${itk-module}Tests
  predicates_test.cxx
)

CreateTestDriver(${itk-module}  "${${itk-module}-Test_LIBRARIES}" "${ITK${itk-module}Tests}")

itk_add_test(NAME IsPointInCircleTest_1
  COMMAND ${itk-module}TestDriver predicates_test
  )
