
foreach(libname BLAS/SRC/libblas.a
    F2CLIBS/libf2c/libf2c.a
    SRC/liblapack.a
    TESTING/MATGEN/libtmglib.a)
  get_filename_component(fname ${libname} NAME)

  configure_file(${clapack_prefix}/${libname}
    ${prereqs}/lib/${fname} COPYONLY)

endforeach(libname)