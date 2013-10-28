file(READ ${AIR_FILE} code)
string(REPLACE "#if defined(_WIN32) || defined(__ECC) /* NrrdIO-hack-002 */"
  "#if defined(_WIN32) || defined(__ECC) || defined(__INTEL_COMPILER) /* NrrdIO-hack-002 */"
  code "${code}" )
file(WRITE ${AIR_FILE} "${code}")
