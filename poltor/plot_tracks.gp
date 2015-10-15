#!/home/palken/usr/bin/gnuplot
#
# Plot single tracks
#
# Adjust 'outdir' below

outdir = 'png'

filedata = 'plp_B.dat'
filestats = 'track_stats.dat.B'

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

l = 0.6
r = 0.0
b = 0.15
vbuffer = 0.3
plotwidth = 4.0
plotheight = 1.0
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

stats filestats;
nplot = STATS_records

load 'xyborder.cfg'
load 'grid.cfg'

do for [idx=0:nplot - 1] {

np = idx + 1

# Retrieve timestamp string
cmd = sprintf('/nfs/satmag_work/palken/corr/common/time2str2 %s', word(tarr, np))
tstr = system(cmd)

# Retrieve longitude of equator crossing
phi = word(phiarr, np) + 0
phistr = sprintf('%.1f', phi)

# Retrieve LT of equator crossing
lt = word(ltarr, np) + 0
lthour = int(lt)
ltmin = int((lt - lthour) * 60)
ltstr = sprintf('%02d:%02d', lthour, ltmin)

# Retrieve mean altitude of track
alt = word(altarr, np) + 0
altstr = sprintf('%.1f km', alt)

cmd = sprintf('/nfs/satmag_work/palken/corr/common/time2str3 %s', word(tarr, np))
fstr = system(cmd)
set out sprintf('%s/plot_%s.png', outdir, fstr)

str = sprintf('Generating plot %d/%d...', np, nplot)
print str

set multiplot layout nrow,ncol

load 'multi_reset.cfg'

set format x ""
set xrange [-60:60]
load 'lines2.cfg'

set ylabel "X residual (nT)"
set title tstr.', {/Symbol \152} = '.phistr.'{/Symbol \260}'.', LT = '.ltstr.', alt = '.altstr.', track '.np

plot filedata us 9:(abs($11) < 50 ? $11 : 1/0) index idx w li lw 4 ti ""

unset title

load 'incrow.cfg'

set ylabel "Z residual (nT)"
plot filedata us 9:(abs($13) < 50 ? $13 : 1/0) index idx w li lw 4 ti ""

load 'incrow.cfg'

set ylabel "F residual (nT)"
plot filedata us 9:14 index idx w li lw 4 ti ""

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set format x "%g"
set format y "%.1t"
set yrange [0:2.5e6]
set label 1 "x10^6" at graph -0.02,1.1

set ylabel "electron density (cm^{-3})"
plot filedata us 9:19 index idx w li lw 4 lt 3 ti ""

unset label 1
unset xlabel
set format x ""
set format y "%g"
set yrange [*:*]

unset multiplot

}
