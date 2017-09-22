/*
 * gridgeo.c
 *
 * Bins (x, y) data into x bins, computes the mean and stddev of
 * the y values in each bin and outputs the results
 *
 * Usage: ./gridgeo [-x xcol] [-y ycol] [-z zcol] [-n numbins] [-N maxdata]
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <common/bin2d.h>

#define MAX_BUFFER  2048

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --xmin | -a xmin          - minimum x\n");
  fprintf(stderr, "\t --xmax | -b xmax          - maximum x\n");
  fprintf(stderr, "\t --ymin | -c ymin          - minimum y\n");
  fprintf(stderr, "\t --ymax | -d ymax          - maximum z\n");
  fprintf(stderr, "\t --xcol | -x xcol          - x column\n");
  fprintf(stderr, "\t --ycol | -y ycol          - y column\n");
  fprintf(stderr, "\t --zcol | -z zcol          - z column\n");
  fprintf(stderr, "\t --xwidth | -u xwidth      - x bin width\n");
  fprintf(stderr, "\t --ywidth | -v ywidth      - y bin width\n");
  fprintf(stderr, "\t --output_file | -o file   - output file\n");
}

int
main(int argc, char *argv[])
{
  char *outfile = "data.grd";
  size_t npts, nx, ny;
  double x, y, z;
  double xwidth, ywidth;
  size_t xcol, ycol, zcol;
  char buffer[MAX_BUFFER];
  char *bufptr;
  double xmin, xmax, ymin, ymax;
  bin2d_workspace *bin2d_p;

  xmin = -180.0;
  xmax = 180.0;
  ymin = -90.0;
  ymax = 90.0;
  xcol = 1;
  ycol = 2;
  zcol = 3;
  xwidth = 5.0;
  ywidth = 5.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "xmin", required_argument, NULL, 'a' },
          { "xmax", required_argument, NULL, 'b' },
          { "ymin", required_argument, NULL, 'c' },
          { "ymax", required_argument, NULL, 'd' },
          { "xcol", required_argument, NULL, 'x' },
          { "ycol", required_argument, NULL, 'y' },
          { "zcol", required_argument, NULL, 'z' },
          { "xwidth", required_argument, NULL, 'u' },
          { "ywidth", required_argument, NULL, 'v' },
          { "output_file", required_argument, NULL, 'o' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:d:o:u:v:x:y:z:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            xmin = atof(optarg);
            break;

          case 'b':
            xmax = atof(optarg);
            break;

          case 'c':
            ymin = atof(optarg);
            break;

          case 'd':
            ymax = atof(optarg);
            break;

          case 'u':
            xwidth = atof(optarg);
            break;

          case 'v':
            ywidth = atof(optarg);
            break;

          case 'x':
            xcol = strtol(optarg, NULL, 0);
            break;

          case 'y':
            ycol = strtol(optarg, NULL, 0);
            break;

          case 'z':
            zcol = strtol(optarg, NULL, 0);
            break;

          case 'o':
            outfile = optarg;
            break;

          case '?':
          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (xcol == ycol)
    {
      fprintf(stderr, "x column cannot equal y column\n");
      exit(1);
    }

  nx = (size_t) round((xmax - xmin) / xwidth);
  ny = (size_t) round((ymax - ymin) / ywidth);

  bin2d_p = bin2d_alloc(xmin, xmax, nx, ymin, ymax, ny);

  fprintf(stderr, "main: x bin size = %g (%zu bins)\n", xwidth, nx);
  fprintf(stderr, "main: y bin size = %g (%zu bins)\n", ywidth, ny);

  fprintf(stderr, "main: reading data points...");

  npts = 0;
  while (fgets(buffer, MAX_BUFFER, stdin) != NULL)
    {
      double val;
      int done = 0;
      int cnt;
      size_t col;

      if (*buffer == '#')
        continue;

      bufptr = buffer;

      col = 1;
      while (!done)
        {
          int n;

          n = sscanf(bufptr, "%lf %n", &val, &cnt);
          if (n <= 0)
            {
              done = 1;
              continue;
            }

          bufptr += cnt;

          if (col == xcol)
            x = val;
          else if (col == ycol)
            y = val;
          else if (col == zcol)
            z = val;

          ++col;
        }

      if (x < xmin || x > xmax || y < ymin || y > ymax)
        continue;

      /* add (x, y, z) to dataset */
      bin2d_add_element(x, y, z, bin2d_p);

      ++npts;
    }

  fprintf(stderr, "done (%zu data points read)\n", npts);

  fprintf(stderr, "main: writing grid to %s...", outfile);
  bin2d_print(outfile, bin2d_p);
  fprintf(stderr, "done\n");

  bin2d_free(bin2d_p);

  return 0;
} /* main() */
