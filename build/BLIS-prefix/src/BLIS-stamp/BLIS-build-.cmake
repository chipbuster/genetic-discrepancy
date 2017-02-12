

set(command "${make}")
execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_FILE "/work/02711/chipbus/star-discrepancy-project/src/build/BLIS-prefix/src/BLIS-stamp/BLIS-build-out.log"
  ERROR_FILE "/work/02711/chipbus/star-discrepancy-project/src/build/BLIS-prefix/src/BLIS-stamp/BLIS-build-err.log"
  )
if(result)
  set(msg "Command failed: ${result}\n")
  foreach(arg IN LISTS command)
    set(msg "${msg} '${arg}'")
  endforeach()
  set(msg "${msg}\nSee also\n  /work/02711/chipbus/star-discrepancy-project/src/build/BLIS-prefix/src/BLIS-stamp/BLIS-build-*.log")
  message(FATAL_ERROR "${msg}")
else()
  set(msg "BLIS build command succeeded.  See also /work/02711/chipbus/star-discrepancy-project/src/build/BLIS-prefix/src/BLIS-stamp/BLIS-build-*.log")
  message(STATUS "${msg}")
endif()
