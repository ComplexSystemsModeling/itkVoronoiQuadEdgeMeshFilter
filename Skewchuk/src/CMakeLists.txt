set(ITKSkewchuk_SRC
predicates.c
)

add_library(${itk-module} ${${itk-module}_SRC})

target_link_libraries(${itk-module}  ${${itk-module}_LIBRARIES})

itk_module_target(${itk-module})

