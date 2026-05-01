
import sys, pickle, netCDF4 as nc, numpy as np
from pathlib import Path
args   = pickle.loads(sys.stdin.buffer.read())
date_s, lat_c, lon_c, half, temps, ri, gi, bi = args
nc_path = Path(temps) / f'date_{date_s}.nc'
with nc.Dataset(str(nc_path)) as ds:
    la = ds['lat'][:]; lo = ds['lon'][:]
    rv = ds['reflectance']; nd = rv.ndim
    ci = int(np.argmin(np.abs(la - lat_c)))
    cj = int(np.argmin(np.abs(lo - lon_c)))
    i0=max(0,ci-half); i1=min(len(la),ci+half)
    j0=max(0,cj-half); j1=min(len(lo),cj+half)
    r=np.array(rv[0,ri,i0:i1,j0:j1] if nd==4 else rv[ri,i0:i1,j0:j1],dtype=float)
    g=np.array(rv[0,gi,i0:i1,j0:j1] if nd==4 else rv[gi,i0:i1,j0:j1],dtype=float)
    b=np.array(rv[0,bi,i0:i1,j0:j1] if nd==4 else rv[bi,i0:i1,j0:j1],dtype=float)
    raw=np.stack([r,g,b],axis=0)
    raw[(raw<-0.01)|(raw>1.2)]=np.nan
    dci=ci-i0; dcj=cj-j0; h=half
    while h>=15:
        s0=max(0,dci-h); s1=min(raw.shape[1],dci+h)
        t0=max(0,dcj-h); t1=min(raw.shape[2],dcj+h)
        if np.isfinite(raw[0,s0:s1,t0:t1]).mean()>=0.30: break
        h=max(15,h*2//3)
    s0=max(0,dci-h); s1=min(raw.shape[1],dci+h)
    t0=max(0,dcj-h); t1=min(raw.shape[2],dcj+h)
    crop=raw[:,s0:s1,t0:t1]
    lat_patch=np.array(la[i0:i1][s0:s1],dtype=float)
    lon_patch=np.array(lo[j0:j1][t0:t1],dtype=float)
    # Flip to N-up (south→north) if stored N→S
    if lat_patch[0] > lat_patch[-1]:
        crop=crop[:,::-1,:]
        lat_patch=lat_patch[::-1]
    def st(a):
        lo2=np.nanpercentile(a,2); hi=np.nanpercentile(a,98)
        return np.clip((a-lo2)/max(hi-lo2,1e-8),0,1)
    rgb=np.stack([st(crop[0]),st(crop[1]),st(crop[2])],axis=-1)
    rgb[~np.isfinite(crop[0])]=0.35
    # star position in the flipped array
    star_row_flipped = (crop.shape[1]-1) - (dci-s0)
    sys.stdout.buffer.write(pickle.dumps((
        rgb,
        float(lat_patch[0]), float(lat_patch[-1]),
        float(lon_patch[0]), float(lon_patch[-1]),
        float(lat_c), float(lon_c),
        float(np.isfinite(crop[0]).mean())
    )))
