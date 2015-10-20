/*
 * log.h
 */

#ifndef INCLUDED_log_h
#define INCLUDED_log_h

#include <stdio.h>
#include <stdarg.h>

#define LOG_MAX_BUFFER    4096

typedef struct
{
  FILE *fp;
  size_t flags;
} log_workspace;

#define LOG_TIMESTAMP    (1 << 0) /* timestamp log entries */
#define LOG_APPEND       (1 << 1)
#define LOG_WRITE        (1 << 2)

/*
 * Prototypes
 */

log_workspace *log_alloc(size_t flags, const char *filename, ...);
void log_free(log_workspace *w);
int log_proc(log_workspace *w, const char *format, ...);

#endif /* INCLUDED_log_h */
