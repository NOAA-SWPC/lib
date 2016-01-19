/*
 * log.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include "log.h"

/*
log_alloc()
  Allocate a log workspace

Inputs: flags    - LOG_xxx (write,append)
        filename - log filename
        ...      - additional arguments if 'filename' contains % symbols
*/

log_workspace *
log_alloc(size_t flags, const char *filename, ...)
{
  log_workspace *w;
  va_list args;
  char filebuf[LOG_MAX_BUFFER + 1];
  int n;
  char fstr[3] = "w";

  w = calloc(1, sizeof(log_workspace));
  if (!w)
    {
      fprintf(stderr, "log_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  va_start(args, filename);

  n = vsnprintf(filebuf, LOG_MAX_BUFFER, filename, args);
  if (n > LOG_MAX_BUFFER)
    {
      fprintf(stderr, "log_alloc: buffer not large enough\n");
      return 0;
    }

  va_end(args);

  if (flags & LOG_APPEND)
    *fstr = 'a';

  w->fp = fopen(filebuf, fstr);
  if (!w->fp)
    fprintf(stderr, "log_alloc: error: %s: %s\n", filebuf, strerror(errno));

  w->flags = flags;

  return w;
} /* log_alloc() */

void
log_free(log_workspace *w)
{
  if (w->fp)
    fclose(w->fp);

  free(w);
} /* log_free() */

/*
log_proc()
  Write a log entry

Return: number of characters written
*/

int
log_proc(log_workspace *w, const char *format, ...)
{
  va_list args;
  int n = 0;

  if (!w)
    return 0;

  va_start(args, format);

  if (w->flags & LOG_TIMESTAMP)
    {
      time_t t = time(0);
      struct tm *tm_p = gmtime(&t);

      n += fprintf(w->fp, "[%04d-%02d-%02d %02d:%02d:%02d] ",
                   tm_p->tm_year + 1900,
                   tm_p->tm_mon + 1,
                   tm_p->tm_mday,
                   tm_p->tm_hour,
                   tm_p->tm_min,
                   tm_p->tm_sec);
    }

  n += vfprintf(w->fp, format, args);

  va_end(args);

  fflush(w->fp);

  return n;
} /* log_proc() */
