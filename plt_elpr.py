# plot electron or proton integrated fluxes.
# >python plt_elpr.py 2016-11-03T00:00 2016-11-06T12:00

import datetime, sys, calendar, requests, io
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import dates

start_time = datetime.datetime.strptime(sys.argv[1], '%Y-%m-%dT%H:%M')
end_time = datetime.datetime.strptime(sys.argv[2], '%Y-%m-%dT%H:%M')
flux = sys.argv[3]
fluxes = ['e2','p1','p2','p3','p4','p5','p6','p7']
if not flux in fluxes:
    print('\
Choose from the followings;\n\
=======================\n\
e2: electron > 2 MeV\n\
p1: proton > 1 MeV\n\
p2: proton > 5 MeV\n\
p3: proton > 10 MeV\n\
p4: proton > 30 MeV\n\
p5: proton > 50 MeV\n\
p6: proton > 60 MeV\n\
p7: proton > 100 MeV\n\
======================\
')
    sys.exit(0)

flux += '_flux_ic'
base_url    = 'https://satdat.ngdc.noaa.gov/sem/goes/data/avg/'
omni_url    = 'https://cdaweb.sci.gsfc.nasa.gov/pub/data/omni/low_res_omni/'
year        = start_time.year
month       = start_time.month
while True:
    print(year,month)
    lastday = calendar.monthrange(year,month)[1]
    filename6 = base_url+\
        '{:04}/{:02}/goes06/csv/g06_eps_5m_{:04}{:02}01_{:04}{:02}{:02}.csv'.format(
        year,month,year,month,year,month,lastday)
    filename7 = base_url+\
        '{:04}/{:02}/goes07/csv/g07_eps_5m_{:04}{:02}01_{:04}{:02}{:02}.csv'.format(
        year,month,year,month,year,month,lastday)
    try:
        r = requests.get(filename7)
        file7 = io.StringIO(r.text.split('data:')[1])
        tmpelpr_df = pd.read_csv(file7, index_col=0, na_values=-99999)
    except:
        r = requests.get(filename6)
        file6 = io.StringIO(r.text.split('data:')[1])
        tmpelpr_df = pd.read_csv(file6, index_col=0, na_values=-99999)

    if (year == start_time.year) & (month == start_time.month):
        elpr_df = tmpelpr_df
    else:
        elpr_df = pd.concat([elpr_df, tmpelpr_df])

    month += 1
    if month > 12:
        month = 1
        year += 1
    if (year == end_time.year) & (month == end_time.month+1):
        break

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
elpr_df.index = pd.to_datetime(elpr_df.index)
targetelpr_df = elpr_df.loc[start_time:end_time]
ax.plot(targetelpr_df.index, targetelpr_df[flux], color='black')
myFmt = dates.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
plt.yscale('log')
ax.set_ylim(1e0,1e7)
ax.set_xlim([start_time,end_time])
plt.xticks(rotation=30, size=15, horizontalalignment='right')
plt.grid()
plt.xlabel('Time (Starting at '+\
        start_time.strftime('%Y-%m-%d %H:%M')+' UT)', size=15)
plt.ylabel('pfu', size=15)
plt.subplots_adjust(left=0.12, right=0.98,bottom=0.1,top=0.95)


