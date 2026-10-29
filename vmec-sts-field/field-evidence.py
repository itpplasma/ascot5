"""Synthetic native VMEC fixture; libascot B_STS vs independent Fourier sums.
Run next to importdata-before.py using the fixed a5py installation.
No measured machine or equilibrium data is used.
"""
import importlib.util
import json
from pathlib import Path
import numpy as np
import netCDF4
import unyt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a5py import Ascot
from a5py.templates.importdata import ImportData
from a5py.testascot.unittests import TestVmecField

fixture = TestVmecField('test_field')
fixture.setUp()
wout = fixture.writewout(-2.0)
with netCDF4.Dataset(wout) as f:
    vals = {k: np.asarray(f[k][:]) for k in ['xm','xn','rmnc','zmns','bsupumnc','bsupvmnc']}
xm, xn = vals['xm'], vals['xn']
ns = vals['rmnc'].shape[0]
surfaces = np.array([.25,.5,.7,.8,.9,.95,.97,.99,1.,1.01,1.02,1.05])
rng = np.random.default_rng(12029)
u, v = rng.uniform(0,2*np.pi,(2,4000))
cos, sin = np.cos(u[:,None]*xm-v[:,None]*xn), np.sin(u[:,None]*xm-v[:,None]*xn)
def radial(key, s, half=False):
    x = s*(ns-1) + (.5 if half else 0)
    i = int(np.clip(np.floor(x),1 if half else 0,ns-2))
    return vals[key][i]*(i+1-x) + vals[key][i+1]*(x-i)
def direct(s):
    rc, zc = radial('rmnc',s), radial('zmns',s)
    r, z = cos@rc, sin@zc
    bu, bv = cos@radial('bsupumnc',s,True), cos@radial('bsupvmnc',s,True)
    b = np.array([bu*(-sin@(xm*rc))+bv*(sin@(xn*rc)),bv*r,
                  bu*(cos@(xm*zc))-bv*(cos@(xn*zc))])
    return r,z,b
spec = importlib.util.spec_from_file_location('old_importdata', 'importdata-before.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
fixed = ImportData.vmec_field
records = []
for label, template in [('Before',old.ImportData.vmec_field),('After',fixed)]:
    ImportData.vmec_field = staticmethod(template)
    fn = Path('synthetic-'+label.lower()+'.h5')
    fn.unlink(missing_ok=True)
    a5 = Ascot(str(fn),create=True)
    a5.data.create_input('vmec_sts',ncfile=wout,phimin=0,phimax=72,nphi=18,ntheta=120,nr=50,nz=50,activate=True)
    d = a5.data.bfield.active.read()
    a5.input_init(bfield=True)
    try:
        for s in surfaces:
            r,z,b = direct(s)
            inside = ((r>d['b_rmin'][0])&(r<d['b_rmax'][0])&(z>d['b_zmin'][0])&(z<d['b_zmax'][0]))
            out = a5.input_eval(r[inside]*unyt.m,v[inside]*unyt.rad,z[inside]*unyt.m,0*unyt.s,'br','bphi','bz','rho')
            got = np.array([o.v for o in out[:3]])
            assert np.isfinite(got).all() and np.isfinite(out[3].v).all()
            err = np.linalg.norm(got-b[:,inside],axis=0)/np.linalg.norm(b[:,inside],axis=0)
            drho = out[3].v-np.sqrt(s)
            records.append(dict(version=label,s=float(s),n=int(inside.sum()),fraction=float(inside.mean()),B_rms=float(np.sqrt(np.mean(err**2))),B_max=float(err.max()),rho_rms=float(np.sqrt(np.mean(drho**2))),rho_max=float(np.max(np.abs(drho)))))
    finally:
        a5.input_free()
ImportData.vmec_field = staticmethod(fixed)
Path('synthetic-field-errors.json').write_text(json.dumps(records,indent=2)+'\n')
fig, axes = plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
for label, color, marker, style in [('Before','#D55E00','o','--'),('After','#0072B2','s','-')]:
    data = [d for d in records if d['version']==label]
    for ax,key in zip(axes,['B_rms','rho_rms']):
        ax.semilogy(surfaces,[d[key] for d in data],label=label,color=color,marker=marker,ls=style,ms=4)
for ax in axes:
    ax.axvline(1,color='0.45',lw=.8,ls=':')
    ax.set_xlabel(r'Normalized toroidal flux $s$ [1]')
    ax.grid(alpha=.2,which='both')
    ax.legend(loc='best',frameon=False)
axes[0].set_ylabel(r'RMS $\|\mathbf{B}_{STS}-\mathbf{B}_{ref}\|/\|\mathbf{B}_{ref}\|$ [1]')
axes[1].set_ylabel(r'RMS $|\rho_{STS}-\sqrt{s}|$ [1]')
fig.suptitle('Synthetic rotating ellipse: ASCOT field evaluation vs direct VMEC sums',fontsize=12)
fig.savefig('before-after-synthetic.png',dpi=180)
print(json.dumps(records,indent=2))
