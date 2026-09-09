if (UNIX)
  # this will allow to use same _DEBUG macro available in both Linux
  # as well as Windows - MSVC environment.
  add_compile_options("$<$<CONFIG:DEBUG>:-D_DEBUG>")
endif (UNIX)
