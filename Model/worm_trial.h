// First shot at a header file, following the example at: 
//  http://stackoverflow.com/questions/15891781/how-to-call-on-a-function-found-on-another-file

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string.h>

// My own additions (to get it to work with sundials 2.6.0)
#include <stdarg.h>     /* va_list, va_start, va_arg, va_end */

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

int worm_trial(const char* pramsFileName);

int write_prams_file(int n, ...);