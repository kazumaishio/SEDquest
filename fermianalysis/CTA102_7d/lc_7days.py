import numpy as np
import pandas as pd
from fermipy.gtanalysis import GTAnalysis
from time import monotonic


#(base)  14:09:41 1_3C454  $ ls result1/
#bexpmap_00.fits     ccube_00.fits       evfile_00.txt       gtbin.par           gtltcube.par        gtsrcmaps.par       ltcube_00.fits
#bexpmap_roi_00.fits config.yaml         ft1_00.fits         gtexpcube2.par      gtselect.par        list.txt            srcmdl_00.xml

name  = 'CTA 102'
sname = 'CTA_102'
binsize = 7

start = monotonic()

# setup the analysis
gta = GTAnalysis('config.yaml', logging={'verbosity': 3})
gta.setup()

# model
gta.free_sources(distance=3.0, pars='norm')

# Free all parameters of isotropic and galactic diffuse components
gta.free_source('galdiff')
gta.free_source('isodiff')

# fit
fit_results = gta.fit()

# save
gta.write_roi(f'fit0_{binsize}days')

# SED
sed = gta.sed(name, use_local_index=True)

'''
# create DataFrame
df2 = pd.DataFrame()
df2['logE_ctr'] = sed['loge_ctr']
df2['flux'] = sed['flux']
df2['spectral_index'] = sed['index']

# save dat
with open(f'result1/{sname}_sed_{b}days.dat', 'w') as f1:
    df2.rename(columns={'logE_ctr': '#logE_ctr'}).to_csv(f1, index=False, float_format='%.12e', sep='\t')

# save csv
with open(f'{sname}_sed_{binsize}days.csv', 'w') as f2:
    df2.to_csv(f2, index=None)
'''

# LC
lc = gta.lightcurve(name, binsz=86400.*binsize, save_bin_data=True, shape_ts_threshold=25.0, use_scaled_srcmap=True, multithread=True, nthread=4, outdir=f'{binsize}days')

print(lc['tmin_mjd'])
print(lc['tmax_mjd'])
print(lc['fit_success'])
print(lc['ts_var'])
print(lc['flux'])
print(lc['eflux'])
print(lc['flux_ul95'])

# create DataFrame
df = pd.DataFrame()
df['tmin_mjd'] = lc['tmin_mjd']
df['tmax_mjd'] = lc['tmax_mjd']
df['flux'] = lc['flux']
df['flux_err'] = lc['flux_err']
df['ts'] = lc['ts']
df['dnde'] = lc['dnde']
df['dnde_err'] = lc['dnde_err']
df['dnde_index'] = lc['dnde_index']

# save dat
with open(f'{sname}_lc_{binsize}days.dat', 'w') as f3:
    df.rename(columns={'tmin_mjd': '#tmin_mjd'}).to_csv(f3, index=False, float_format='%.12e', sep='\t')

# save csv
with open(f'{sname}_lc_{binsize}days.csv', 'w') as f4:
    df.to_csv(f4, index=None)


end = monotonic()
print(f"Czas kompilacji: {(end-start)/60./60.} h")
