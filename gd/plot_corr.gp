#!/home/palken/usr/bin/gnuplot
#
# Plot tracks
#
# Adjust 'outdir' below

outdir = 'plots'

filepeak = 'peak.dat'
filecorr = 'corr.dat'

# Read timestamps into array
set term dumb
tarr=""
phiarr=""
ltarr=""
neqd1arr=""
ne1arr=""
neqd2arr=""
ne2arr=""
Fqd1arr=""
F1arr=""
Fqd2arr=""
F2arr=""
corr1arr=""
corr2arr=""
plot filecorr us ( tarr=tarr.stringcolumn(1).' ', \
                   phiarr=phiarr.stringcolumn(3).' ', \
                   ltarr=ltarr.stringcolumn(4).' ', \
                   neqd1arr=neqd1arr.stringcolumn(7).' ', \
                   ne1arr=ne1arr.stringcolumn(8).' ', \
                   neqd2arr=neqd2arr.stringcolumn(9).' ', \
                   ne2arr=ne2arr.stringcolumn(10).' ', \
                   Fqd1arr=Fqd1arr.stringcolumn(11).' ', \
                   F1arr=F1arr.stringcolumn(12).' ', \
                   Fqd2arr=Fqd2arr.stringcolumn(13).' ', \
                   F2arr=F2arr.stringcolumn(14).' ', \
                   corr1arr=corr1arr.stringcolumn(15).' ', \
                   corr2arr=corr2arr.stringcolumn(16).' ', $1):2

nrow = 3
ncol = 1

load 'multi_default.cfg'

l = 1.0
r = 0.4
b = 0.2
vbuffer = 0.3
plotwidth = 5.0
plotheight = 1.2
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

stats filecorr;
nplot = STATS_records

load 'grid.cfg'
load 'lines2.cfg'

set xrange [-60:60]

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

cmd = sprintf('/nfs/satmag_work/palken/corr/common/time2str3 %s', word(tarr, np))
fstr = system(cmd)
set out sprintf('%s/plot_%s.png', outdir, fstr)

str = sprintf('Generating plot %d/%d...', np, nplot)
print str

set multiplot layout nrow,ncol

load 'multi_reset.cfg'

set ylabel "electron density (cm^{-3})"

set title tstr.', {/Symbol \152} = '.phistr.'{/Symbol \260}'.', LT = '.ltstr.', track '.np
load 'xy2border.cfg'
set format x ""
set yrange [0:2.5e6]
set y2range [-4e4:4e4]

neqd1 = word(neqd1arr, np) + 0
ne1 = word(ne1arr, np) + 0
neqd2 = word(neqd2arr, np) + 0
ne2 = word(ne2arr, np) + 0

plot filepeak us 1:2 index idx w li lw 4 ti "N_e", \
     filepeak us 1:4 index idx w li lw 2 axes x1y2 ti "dN_e", \
     filepeak us 1:5 index idx w li lw 2 axes x1y2 ti "smoothed dN_e", \
     '+' us (neqd1):(ne1) w p ps 2 pt 6 lc rgb "red" ti "", \
     '+' us (neqd2):(ne2) w p ps 2 pt 6 lc rgb "red" ti ""

set format x ""
unset title

load 'incrow.cfg'

set yrange [*:*]
set y2range [-0.1:0.1]

Fqd1 = word(Fqd1arr, np) + 0
F1 = word(F1arr, np) + 0
Fqd2 = word(Fqd2arr, np) + 0
F2 = word(F2arr, np) + 0

plot filepeak us 1:3 index idx w li lw 4 ti "F", \
     filepeak us 1:6 index idx w li lw 2 axes x1y2 ti "dF", \
     filepeak us 1:7 index idx w li lw 2 axes x1y2 ti "smoothed dF", \
     '+' us (Fqd1):(F1) w p ps 2 pt 6 lc rgb "red" ti "", \
     '+' us (Fqd2):(F2) w p ps 2 pt 6 lc rgb "red" ti ""

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"
set format x "%g"

set yrange [0:2.5e6]
set y2range [*:*]

corr1 = word(corr1arr, np) + 0
corr2 = word(corr2arr, np) + 0
corr1str = sprintf('north peak corr = %.2f', corr1)
corr2str = sprintf('south peak corr = %.2f', corr2)

set label 1 corr1str at graph 0.01,0.9
set label 2 corr2str at graph 0.01,0.8

plot filepeak us 1:2 index idx w li lw 4 ti "N_e", \
     filepeak us 1:3 index idx w li lw 4 axes x1y2 ti "F", \
     filepeak us 1:8 index idx w li lt 4 lw 4 ti "fit data 1", \
     filepeak us 1:9 index idx w li lt 5 lw 4 ti "fit data 2", \
     filepeak us 1:10 index idx w li lt 4 lw 4 axes x1y2 ti "", \
     filepeak us 1:11 index idx w li lt 5 lw 4 axes x1y2 ti "", \
     '+' us (neqd1):(ne1) w p ps 2 pt 6 lc rgb "red" ti "", \
     '+' us (neqd2):(ne2) w p ps 2 pt 6 lc rgb "red" ti "", \
     '+' us (Fqd1):(F1) w p axes x1y2 ps 2 pt 6 lc rgb "red" ti "", \
     '+' us (Fqd2):(F2) w p axes x1y2 ps 2 pt 6 lc rgb "red" ti ""

unset label 1
unset label 2
unset xlabel

unset multiplot

}
