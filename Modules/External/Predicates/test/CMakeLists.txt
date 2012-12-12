itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(ITK${itk-module}Tests
  IsPointInCircle.cxx
)

CreateTestDriver(${itk-module}  "${${itk-module}-Test_LIBRARIES}" "${ITK${itk-module}Tests}")



# in test one we take one point of the square and move it AWAY from the center
# of the circle. The point is outside, the test fail as expected.
itk_add_test(NAME TestIsPointInCircle_1
  COMMAND ${itk-module}TestDriver IsPointInCircle -0.2 -0.2 0
  )
set_tests_properties(
    TestIsPointInCircle_1 
    PROPERTIES WILL_FAIL "TRUE"
    )

# in test two we take one point of the square and move it TOWARD the center
# of the circle. The point is then inside, the test passes.
itk_add_test(NAME TestIsPointInCircle_2
  COMMAND ${itk-module}TestDriver IsPointInCircle 0.2 0.2 0
  )

# in test 3 we move the point only sightly, the point is now outside, but
# rounding errors make the test pass. This is the example of the limitation 
# of float arithmetic
itk_add_test(NAME TestIsPointInCircle_3
  COMMAND ${itk-module}TestDriver IsPointInCircle 0 -0.00000000000000002 0
  )

# in test 4, we move the point as in test 3 but we use exact geometrical
# predicate. This time, the test as the correct behavior.
itk_add_test(NAME TestIsPointInCircle_4
  COMMAND ${itk-module}TestDriver IsPointInCircle 0 -0.00000000000000002 1
  )
set_tests_properties(
    TestIsPointInCircle_4 
    PROPERTIES WILL_FAIL "TRUE"
    )
