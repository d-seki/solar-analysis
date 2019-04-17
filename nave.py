import numpy as np
from sunpy.map import Map
from scipy import ndimage, optimize, interpolate
import matplotlib.pyplot as plt


"""
Nonlinear Affine Velocity Estimator (NAVE) based on Chae & Sakurai 2004
Differential Affine Velocity Estimator (DAVE) based on Chae & Sakurai 2004
Local Correlation Tracking based on Chae & Sakurai 2004

Keywords------------
prvs:   before image
nxt:    after image
pos:    position where velocity is calculated
        pos[0]:x
        pos[1]:y
lfs:    window size = 2*lfs+1
        lfs = 1: windw size 3x3
        lfs = 2: windw size 5x5

"""


class nave:
    def __init__(self, prvs, nxt, lfs, mode='lct'):
        self.lfs            = lfs
        self.mode           = mode
        self.prvsshape      = prvs.shape
        self.nxtshape       = nxt.shape
        self.Lx             = prvs.shape[1]
        self.Ly             = prvs.shape[0]
        x_array             = np.arange(self.Lx)
        y_array             = np.arange(self.Ly)
        x_grid, y_grid      = np.meshgrid(x_array, y_array)
        x_grid              = x_grid.reshape(self.Lx, self.Ly,1,1)
        y_grid              = y_grid.reshape(self.Lx, self.Ly,1,1)
        self.xy_array       = np.concatenate((x_grid, y_grid), axis=-2)
        self.rltvxy_array   = self.xy_array
        self.wndw           = np.zeros((self.Lx, self.Ly))
        self.prmnum         = {'nave':6,'dave':6,'lct':2}
        self.Bfr            = interpolate.interp2d(
                                x_array, y_array, prvs,
                                kind='cubic', copy=True, bounds_error=False,
                                fill_value=0)
        self.Aftr           = interpolate.interp2d(
                                x_array, y_array, nxt,
                                kind='cubic', copy=True, bounds_error=False,
                                fill_value=0)

    def residual_nave(self, UV):
        u0 = np.array([
            [UV[0]],
            [UV[1]]])
        W = np.array([
            [UV[2], UV[3]],
            [UV[4], UV[5]]])
        nu = -(UV[2]+UV[5])
        W_dot_rltvxy_array = np.tensordot(self.rltvxy_array,W, axes=[2,1])
        W_dot_rltvxy_array = W_dot_rltvxy_array.reshape((self.Lx,self.Ly,2,1))
        RtrjctryBfr = (self.xy_array-0.5*u0-0.5*W_dot_rltvxy_array)[
                        self.wndw==1,:,:]
        RtrjctryAftr = (self.xy_array+0.5*u0+0.5*W_dot_rltvxy_array)[
                        self.wndw==1,:,:]
        rsdl_array =\
            np.exp(-0.5*nu)*self.Aftr(RtrjctryAftr[:,0,0],RtrjctryAftr[:,1,0])-\
            np.exp(0.5*nu)*self.Bfr(RtrjctryBfr[:,0,0],RtrjctryBfr[:,1,0])
        return (rsdl_array*rsdl_array).sum()

    def residual_dave(self, UV):
        # NOTICE: just the case in which nu is small enough
        u0 = np.array([
            [UV[0]],
            [UV[1]]])
        W = np.array([
            [UV[2], UV[3]],
            [UV[4], UV[5]]])
        nu = -(UV[2]+UV[5])
        W_dot_rltvxy_array = np.tensordot(self.rltvxy_array,W, axes=[2,1])
        W_dot_rltvxy_array = W_dot_rltvxy_array.reshape((self.Lx,self.Ly,2,1))
        RtrjctryBfr = (self.xy_array-0.5*u0-0.5*W_dot_rltvxy_array)[
                        self.wndw==1,:,:]
        RtrjctryAftr = (self.xy_array+0.5*u0+0.5*W_dot_rltvxy_array)[
                        self.wndw==1,:,:]
        rsdl_array =\
            (1-0.5*nu)*self.Aftr(RtrjctryAftr[:,0,0],RtrjctryAftr[:,1,0])-\
            (1-0.5*nu)*self.Bfr(RtrjctryBfr[:,0,0],RtrjctryBfr[:,1,0])
        return (rsdl_array*rsdl_array).sum()

    def residual_lct(self, UV):
        u0 = np.array([
            [UV[0]],
            [UV[1]]])
        RtrjctryBfr = (self.xy_array - 0.5*u0)[
                        self.wndw==1,:,:]
        RtrjctryAftr = (self.xy_array + 0.5*u0)[
                        self.wndw==1,:,:]
        rsdl_array =\
            self.Aftr(RtrjctryAftr[:,0,0], RtrjctryAftr[:,1,0])-\
            self.Bfr(RtrjctryBfr[:,0,0], RtrjctryBfr[:,1,0])
        return (rsdl_array*rsdl_array).sum()

    def compute(self):
        if self.prvsshape != self.nxtshape:
            print('Shapes of 2 images must be equivalent.')
            return
        if self.mode == 'nave':
            residual_func = self.residual_nave
        elif self.mode == 'dave':
            residual_func = self.residual_dave
        elif self.mode == 'lct':
            residual_func = self.residual_lct
        else:
            print('Mode must be nave, dave, or lct.')
            return
        params = np.zeros((self.Lx, self.Ly, self.prmnum[self.mode]))
        intl_prms  = [0.1,0.1,0.1,0.1,0.1,0.1][:self.prmnum[self.mode]]
        for x0 in range(self.lfs, self.Lx-self.lfs, 2*self.lfs+1):
            print(x0,self.Lx-self.lfs)
            for y0 in range(self.lfs, self.Ly-self.lfs, 2*self.lfs+1):
                xy0 = np.array([
                        [x0],
                        [y0]])
                self.rltvxy_array = self.xy_array - xy0
                self.wndw = np.zeros((self.Lx, self.Ly))
                self.wndw[int(y0-self.lfs):1+int(y0+self.lfs),
                            int(x0-self.lfs):1+int(x0+self.lfs)] = 1.
                sol = optimize.minimize(residual_func, intl_prms)
                # NOTICE: following expression can be applied only for LCT.
                params[y0-self.lfs:y0+self.lfs+1,
                        x0-self.lfs:x0+self.lfs+1,
                        :] = sol.x
        return params

