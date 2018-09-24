#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "Exception.h"


#ifdef __GNUC__
#include <execinfo.h>
// Content from glibc manual:
/* Obtain a backtrace and print it to out. */
void printTrace2(FILE *out, const char *file, int line) {
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);

  fprintf(out, "Call stack from %s:%d:\n", file, line);

  for (i = 0; i < size; i++)
    fprintf(out, "    %s\n", strings[i]);

  free(strings);
}

// Content from http://tombarta.wordpress.com:
#include <string.h>
#include <cxxabi.h>
void printTrace(FILE *out, const char *file, int line) {
  const size_t max_depth = 100;
  size_t stack_depth;
  void *stack_addrs[max_depth];
  char **stack_strings;

  stack_depth = backtrace(stack_addrs, max_depth);
  stack_strings = backtrace_symbols(stack_addrs, stack_depth);

  fprintf(out, "Call stack from %s:%d:\n", file, line);

  for (size_t i = 1; i < stack_depth; i++) {
      size_t sz = 200; // just a guess, template names will go much wider
      char *function = (char *)malloc(sz);
      char *begin = 0, *end = 0;
      // find the parentheses and address offset surrounding the mangled name
      for (char *j = stack_strings[i]; *j; ++j) {
        if (*j == '(') {
          begin = j;
        }
        else if (*j == '+') {
          end = j;
        }
      }
      if (begin && end) {
        *begin++ = 0;
        *end = 0;
        // found our mangled name, now in [begin, end)

        int status;
        char *ret = abi::__cxa_demangle(begin, function, &sz, &status);
        if (ret) {
          // return value may be a realloc() of the input
          function = ret;
        } else {
          // demangling failed, just pretend it's a C function with no args
          strncpy(function, begin, sz);
          strncat(function, "()", sz);
          function[sz-1] = 0;
        }
        fprintf(out, "    %s:%s\n", stack_strings[i], function);
      } else {
          // didn't find the mangled name, just print the whole line
          fprintf(out, "    %s\n", stack_strings[i]);
      }
      free(function);
  }
  free(stack_strings); // malloc()ed by backtrace_symbols
  fflush(out);
}
#endif


// code adopted from printf manpage
Exception::Exception(const char *fmt, ...) throw() {
#ifdef __GNUC__
  printTrace(stderr, "MyException.cpp", 87);
#endif

  _buffer = NULL;

  int n, size = 256;
  va_list ap;

  while (1) {
    if ((_buffer = new char[size]) == NULL)
      return;

    /* Try to print in the allocated space. */
    va_start(ap, fmt);
    n = vsnprintf(_buffer, size, fmt, ap);
    va_end(ap);

    /* If that worked, return the string. */
    if(n > -1 && n < size)
      return;

    /* Else try again with more space. */
    if(n > -1)    /* glibc 2.1 */
        size = n+1; /* precisely what is needed */
    else           /* glibc 2.0 */
        size *= 2;  /* twice the old size */

    delete[] _buffer;
  }
}

Exception::~Exception() throw() {
  if(_buffer) {
    delete[] _buffer;
  }
}

const char* Exception::what() const throw() {
  if(_buffer) {
    return _buffer;
  } else {
    return "Unknown exception!\n";
  }
}
