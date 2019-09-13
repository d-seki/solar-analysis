import datetime, sys, calendar, requests, io
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import dates

start_time = datetime.datetime.strptime(sys.argv[1], '%Y-%m-%dT%H:%M')
end_time = datetime.datetime.strptime(sys.argv[2], '%Y-%m-%dT%H:%M')
flux = sys.argv[3]
fluxes = ['Kp','Dst']
if not flux in fluxes:
    print('\
Choose from the followings;\n\
=======================\n\
Kp: Planetary K-index\n\
Dst: Dst index\n\
======================\
')
    sys.exit(0)

year        = start_time.year
while True:
    print(year)
    tmpOmni_df = pd.read_csv(
        './omni/omni2_{:04}.csv'.format(year),header=None,
        na_values=[99,999,999.9,9999,99999.99,9999999,9.999,99.99,999.99,99999])
    timeindex = pd.DatetimeIndex(
            start=datetime.datetime(year,1,1,0),
            freq='H',
            end=datetime.datetime(year,12,31,23))
    tmpF107VswBKpDst_df = pd.DataFrame({
            'F107'   :tmpOmni_df[50],
            'Vsw'    :tmpOmni_df[24],
            'B'      :tmpOmni_df[9],
            'Bx'     :tmpOmni_df[12],
            'By'     :tmpOmni_df[13],
            'Bz'     :-tmpOmni_df[14],
            'Kp'     :tmpOmni_df[38],
            'Dst'    :-tmpOmni_df[40]})
    tmpF107VswBKpDst_df.index = timeindex

    if (year == start_time.year):
        kpdst_df = tmpF107VswBKpDst_df
    else:
        kpdst_df = pd.concat([kpdst_df, tmpF107VswBKpDst_df])

    year += 1
    if year > end_time.year:
        break

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)
targetkpdst_df = kpdst_df.loc[start_time:end_time]
if flux == 'Dst':
    ax.plot(targetkpdst_df.index, -targetkpdst_df[flux], color='black')
myFmt = dates.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
ax.set_ylim(-300,0)
ax.set_xlim([start_time,end_time])
plt.xticks(rotation=30, size=15, horizontalalignment='right')
plt.grid()
plt.xlabel('Time (Starting at '+\
        start_time.strftime('%Y-%m-%d %H:%M')+' UT)', size=15)
plt.ylabel('Dst [nT]', size=15)
plt.subplots_adjust(left=0.12, right=0.98,bottom=0.1,top=0.95)


