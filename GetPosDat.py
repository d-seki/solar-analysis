import numpy as np
import astropy.units as u
from box import plot_box
import sys, glob, datetime
from sunpy.map import Map
import matplotlib.pyplot as plt
from prep import smartprep, url_aiaprep


class GetPosDat:
    def __init__(self, img, color, mark):
        self.img = img
        self.c = color
        self.mark = mark
        self.posx = []
        self.posy = []
        self.cid = img.figure.canvas.mpl_connect('button_press_event', self)
        self.bl = [img.axes.get_xlim()[0], img.axes.get_ylim()[0]]
        self.tr = [img.axes.get_xlim()[1], img.axes.get_ylim()[1]]
    def __call__(self, event):
        print('you pressed', event)
        if event.inaxes!=self.img.axes: return
        self.posx.append(event.xdata)
        self.posy.append(event.ydata)
        self.img.axes.plot(self.posx, self.posy, c=self.c, marker=self.mark)
    def end(self):
        self.img.figure.canvas.mpl_disconnect(self.cid)
    def length(self):
        dposx = np.diff(self.posx)
        dposy = np.diff(self.posy)
        return np.sqrt(dposx**2 + dposy**2).sum()

#fig = plt.figure(figsize=(7,7))
#img = fitsmap.plot()
#fig.subplots_adjust(left=0.2,right=0.95)
# crop the image
#getposdat = GetPosDat(img)
# take the positions
#getposdat.end()
#fig.savefig('./webpage/img/'+utime.strftime('%y%m%d%H')+'_sub.pdf')
#fig.clear()
#img = fitsmap.plot()
#plot_box(img,getposdat.bl, getposdat.tr, c='white')
#fig.savefig('./webpage/img/'+utime.strftime('%y%m%d%H')+'_full.pdf')
#print((((getposdat.length()*u.pix)*fitsmap.scale[0]).to(u.rad)*fitsmap.coordinate_frame.observer.radius.to(u.Mm)).value)


