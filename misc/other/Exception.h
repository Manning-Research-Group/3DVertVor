#ifndef MY_EXCEPTION_H
#define MY_EXCEPTION_H

#include <stdlib.h>
#include <exception>

#ifdef __GNUC__
#include <stdio.h>
extern void printTrace(FILE *out, const char *file, int line);
#endif

class Exception : public std::exception {
public:
  Exception(const char *fmt=NULL, ...) throw();
  virtual ~Exception() throw();
  virtual const char* what() const throw();
private:
  char *_buffer;
};

#endif

