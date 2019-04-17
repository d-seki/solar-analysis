# plot GOES soft X-ray flux
# >python plt_goes.py 2016-11-03T00:00 2016-11-06T12:00

import datetime, sys
from sunpy.time import TimeRange
from sunpy.net import  Fido, attrs
from sunpy.timeseries import TimeSeries
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import dates


# enter the start and end times as YYYY-mm-ddTHH:MM
start_time = datetime.datetime.strptime(sys.argv[1], '%Y-%m-%dT%H:%M')
end_time = datetime.datetime.strptime(sys.argv[2], '%Y-%m-%dT%H:%M')
tr = TimeRange([start_time, end_time])
results = Fido.search(attrs.Time(tr), attrs.Instrument('XRS'))
files = Fido.fetch(results)
files.sort()
goes = TimeSeries(files)

for i in range(len(goes)):
    if i == 0:
        goesdf = goes[i].data
    else:
        goesdf = pd.concat([goesdf, goes[i].data])

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)
ax.plot(goesdf.index, goesdf.xrsb, color='black')
ax.plot(goesdf.index, goesdf.xrsa, color='gray')
plt.yscale('log')
plt.ylim([0,1e-2])
plt.xticks(rotation=30, size=15, horizontalalignment='right')
plt.yticks(size=15)
plt.xlim([start_time,end_time])
myFmt = dates.DateFormatter('%H:%M')
plt.gca().xaxis.set_major_formatter(myFmt)
plt.title('GOES Xray Flux', size=20)
plt.xlabel('Time (Starting at '+\
        start_time.strftime('%Y-%m-%d %H:%M')+' UT)', size=15)
plt.ylabel('Watts $m^{-2}$', size=15)
plt.grid()
plt.subplots_adjust(left=0.12, right=0.98,bottom=0.22)
