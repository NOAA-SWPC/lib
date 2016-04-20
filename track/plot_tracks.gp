#!/home/palken/usr/bin/gnuplot
#
# Plot tracks
#
# Adjust 'outdir' below

outdir = 'plots'

filedata = 'track_data.dat'
filestats = 'track_stats.dat'

# Read timestamps into array
set term dumb
tarr=""
phiarr=""
ltarr=""
altarr=""
plot filestats us ( tarr=tarr.stringcolumn(1).' ', \
                    phiarr=phiarr.stringcolumn(2).' ', \
                    ltarr=ltarr.stringcolumn(3).' ', \
                    altarr=altarr.stringcolumn(4).' ', $1):2

nrow = 4
ncol = 1

load 'multi_default.cfg'

l = 0.7
r = 0.4
b = 0.2
vbuffer = 0.3
plotwidth = 5.0
plotheight = 1.2
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

stats filestats;
nplot = STATS_records

load 'grid.cfg'
load 'lines2.cfg'
load 'xy2bordercol.cfg'
set format x ""

do for [idx=0:nplot - 1] {

np = idx + 1

# Retrieve timestamp string
cmd = sprintf('/data/palken/lib/common/time2str2 %s', word(tarr, np))
tstr = system(cmd)

# Retrieve longitude of equator crossing
phi = word(phiarr, np) + 0
phistr = sprintf('%.1f', phi)

# Retrieve LT of equator crossing
lt = word(ltarr, np) + 0
lthour = int(lt)
ltmin = int((lt - lthour) * 60)
ltstr = sprintf('%02d:%02d', lthour, ltmin)

# Retrieve altitude of track
alt = word(altarr, np) + 0
altstr = sprintf('altitude = %.1f km', alt)

cmd = sprintf('/data/palken/lib/common/time2str3 %s', word(tarr, np))
fstr = system(cmd)
set out sprintf('%s/plot_%s.png', outdir, fstr)

str = sprintf('Generating plot %d/%d...', np, nplot)
print str

set multiplot layout nrow,ncol

load 'multi_reset.cfg'

set xrange [-50:50]
set xtics 10
set y2range [0:2e6]
set y2label "electron density (cm^{-3})"

set title tstr.', {/Symbol \152} = '.phistr.'{/Symbol \260}'.', LT = '.ltstr.', '.altstr.', track '.np

set ylabel "F residual (nT)"
plot filedata us 8:12 index idx w li lw 4 ti "", \
     filedata us 8:13 index idx w li lw 4 axes x1y2 ti ""

unset title

load 'incrow.cfg'

set ylabel "X residual (nT)"
plot filedata us 8:9 index idx w li lw 4 ti "", \
     filedata us 8:13 index idx w li lw 4 axes x1y2 ti ""

load 'incrow.cfg'

set ylabel "Y residual (nT)"
plot filedata us 8:10 index idx w li lw 4 ti "", \
     filedata us 8:13 index idx w li lw 4 axes x1y2 ti ""

load 'incrow.cfg'

set format x "%g"
set xlabel "QD latitude (degrees)"
set ylabel "Z residual (nT)"

plot filedata us 8:11 index idx w li lw 4 ti "", \
     filedata us 8:13 index idx w li lw 4 axes x1y2 ti ""

set format x ""
unset xlabel

unset multiplot

}
