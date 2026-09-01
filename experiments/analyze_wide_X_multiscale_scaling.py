#!/usr/bin/env python3
"""
Wide-X scaling test for the local multiscale Haar representation.

This script samples prime-square geometry near X=2e7,...,2e12,
fits the short-interval covariance kernel with local Haar scales, tests
coarse L/H scaling and logarithmic rho1 amplitude, and records a negative
fine-spectrum transfer test. It does not count actual primes up to 1e12.
"""
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.sparse import coo_matrix
from scipy.stats import linregress

ROOT=Path(__file__).resolve().parent
X_TARGETS=[2e7,2e8,2e9,2e10,2e11,2e12]
PANEL_M=128
OFFSETS=[-128,0,128]
STEPS_PER_OCTAVE=16
CELLS_PER_MEDIAN_H=64
MAX_LAG=20

def primes_upto(n):
    z=np.ones(n+1,dtype=bool); z[:2]=False
    for p in range(2,math.isqrt(n)+1):
        if z[p]: z[p*p:n+1:p]=False
    return np.flatnonzero(z)

def Delta(t):
    t=float(t)
    return .5*((t+1)*math.log(t+1)-2*t*math.log(t)+(0.0 if t==1 else (t-1)*math.log(t-1)))

DELTA=np.array([Delta(k) for k in range(1,MAX_LAG+1)])
DELTA_NORM=DELTA/DELTA[0]

def V(h,X):
    return 0.0 if h<=0 else h*math.log(X/h)

def pair_rho(A,B,i,j):
    a,b=float(A[i]),float(B[i]); c,d=float(A[j]),float(B[j])
    centers=.5*(A+B)
    Xloc=math.sqrt(float(centers[i])*float(centers[j]))
    cov=.5*(V(d-a,Xloc)+V(c-b,Xloc)-V(c-a,Xloc)-V(d-b,Xloc))
    return cov/math.sqrt(V(b-a,Xloc)*V(d-c,Xloc))

def log_scales(ncell,steps):
    vals=[]; x=1.0; xmax=math.log2(max(2,ncell))
    while x<=xmax+1e-9:
        vals.append(int(round(2**x))); x+=1/steps
    return np.array(sorted(set(v for v in vals if 2<=v<=ncell)),dtype=int)

def haar_basis_adaptive(A,B):
    Hmed=float(np.median(B-A)); cell=max(1,int(round(Hmed/CELLS_PER_MEDIAN_H)))
    origin=int(A[0]); starts=(A-origin).astype(np.int64); ends=(B-origin).astype(np.int64)
    N=int(ends[-1]); M=len(A)
    regular=np.arange(0,N+cell,cell,dtype=np.int64); regular=regular[regular<=N]
    if len(regular)==0 or regular[-1]!=N: regular=np.append(regular,N)
    bounds=np.unique(np.concatenate([regular,starts,ends,[0,N]])).astype(np.int64)
    ss=bounds[:-1]; ee=bounds[1:]; sl=ee-ss
    si=np.searchsorted(ends,ss,side='right').astype(np.int32)
    sc=(ss//cell).astype(np.int64); ilen=(ends-starts).astype(float)
    ncell=int(math.ceil(N/cell)); scales=log_scales(ncell,STEPS_PER_OCTAVE)
    rows=[]
    for L in scales:
        shifts=sorted(set([0,int(round(L/4)),int(round(L/2)),int(round(3*L/4))]))
        cov=np.zeros(MAX_LAG+1)
        for shift in shifts:
            u=sc+shift; block=u//L; pos=u%L; sign=np.where(pos<L/2,1.0,-1.0)
            data=(sl.astype(float)/ilen[si])*sign
            W=coo_matrix((data,(si,block)),shape=(M,int(block.max())+1)).tocsr()
            cov[0]+=float(np.mean(np.asarray(W.multiply(W).sum(axis=1)).ravel()))
            for k in range(1,MAX_LAG+1):
                cov[k]+=float(np.mean(np.asarray(W[:-k].multiply(W[k:]).sum(axis=1)).ravel()))
        rows.append(cov/len(shifts))
    return cell,scales,np.asarray(rows).T

def fit_panel(A,B):
    cell,scales,basis=haar_basis_adaptive(A,B)
    theory=np.array([np.mean([pair_rho(A,B,i,i+k) for i in range(len(A)-k)]) for k in range(1,MAX_LAG+1)])
    target=np.r_[1.0,theory]
    lw=np.ones(MAX_LAG+1); lw[0]=3.0; lw[1:6]=2.0
    v,res=nnls(basis*lw[:,None],target*lw)
    fit=basis@v; corr=fit[1:]/fit[0]; contribution=v*basis[0]
    Hmed=float(np.median(B-A)); Xmed=float(np.median((A+B)/2)); ratios=scales*cell/Hmed
    total=float(contribution.sum())
    def band(lo,hi):
        m=(ratios>=lo)&(ratios<hi)
        return float(contribution[m].sum()/total)
    return {'cell':cell,'scales':scales,'basis':basis,'v':v,'contribution':contribution,
            'theory':theory,'fit_corr':corr,'fit_residual':float(res),
            'kernel_rmse':float(np.sqrt(np.mean((corr-theory)**2))),
            'shape_rmse_vs_Delta':float(np.sqrt(np.mean((corr/corr[0]-DELTA_NORM)**2))),
            'median_H':Hmed,'median_X':Xmed,'X_over_H':Xmed/Hmed,
            'L_over_H_lt_2':band(0,2),'L_over_H_2_4':band(2,4),
            'L_over_H_4_8':band(4,8),'L_over_H_ge_8':band(8,np.inf)}

P_MAX=int(math.sqrt(max(X_TARGETS))*1.05)+5000
primes=primes_upto(P_MAX)
rows=[]; details={}
for Xt in X_TARGETS:
    i0=int(np.argmin(np.abs(primes-math.sqrt(Xt))))
    for off in OFFSETS:
        j=i0+off-PANEL_M//2
        ps=primes[j:j+PANEL_M+1].astype(np.int64)
        A=ps[:-1]**2; B=ps[1:]**2; f=fit_panel(A,B)
        panel_id=f'X{Xt:.0e}_off{off:+d}'; details[panel_id]=(A,B,f)
        rows.append({'panel_id':panel_id,'target_X':Xt,'offset_primes':off,
                     'p_center':int(ps[PANEL_M//2]),'median_X':f['median_X'],
                     'median_H':f['median_H'],'X_over_H':f['X_over_H'],'cell':f['cell'],
                     'kernel_rmse':f['kernel_rmse'],'shape_rmse_vs_Delta':f['shape_rmse_vs_Delta'],
                     'L_over_H_lt_2':f['L_over_H_lt_2'],'L_over_H_2_4':f['L_over_H_2_4'],
                     'L_over_H_4_8':f['L_over_H_4_8'],'L_over_H_ge_8':f['L_over_H_ge_8']})

panels=pd.DataFrame(rows)
panels['L_over_H_lt_4']=panels['L_over_H_lt_2']+panels['L_over_H_2_4']
panels.to_csv(ROOT/'wide_X_multiscale_scaling_panels.csv',index=False,encoding='utf-8-sig')

amp=[]
for panel_id,(A,B,f) in details.items():
    scaled=f['theory'][0]*math.log(f['X_over_H'])
    amp.append({'panel_id':panel_id,'median_X':f['median_X'],'median_H':f['median_H'],
                'X_over_H':f['X_over_H'],'rho1_theory':float(f['theory'][0]),
                'rho1_times_log_X_over_H':float(scaled),
                'difference_from_minus_log2':float(scaled+math.log(2))})
amplitude=pd.DataFrame(amp)
amplitude.to_csv(ROOT/'wide_X_rho1_log_scaling.csv',index=False,encoding='utf-8-sig')

grouped=panels.groupby('target_X').agg({'median_X':['mean','std'],'median_H':['mean','std'],
    'X_over_H':['mean','std'],'L_over_H_lt_4':['mean','std'],'L_over_H_4_8':['mean','std'],
    'L_over_H_ge_8':['mean','std'],'kernel_rmse':['mean','std']}).reset_index()
grouped.columns=['_'.join([str(x) for x in c if str(x)]) for c in grouped.columns.to_flat_index()]
grouped.to_csv(ROOT/'wide_X_multiscale_scaling_grouped.csv',index=False,encoding='utf-8-sig')

logX=np.log10(panels['median_X'].to_numpy(float)); reg=[]
for band_name in ['L_over_H_lt_4','L_over_H_4_8','L_over_H_ge_8']:
    r=linregress(logX,panels[band_name])
    reg.append({'band':band_name,'mean_fraction':float(panels[band_name].mean()),
                'sd_fraction':float(panels[band_name].std(ddof=1)),
                'min_fraction':float(panels[band_name].min()),'max_fraction':float(panels[band_name].max()),
                'slope_per_decade':float(r.slope),'slope_stderr':float(r.stderr),
                'p_value_slope_zero':float(r.pvalue),'R2':float(r.rvalue**2)})
regression=pd.DataFrame(reg)
regression.to_csv(ROOT/'wide_X_multiscale_scaling_regression.csv',index=False,encoding='utf-8-sig')

# Fine universal-spectrum transfer diagnostic.
R_REF=200.0
H=int(round(20_000_000/R_REF)); a0=int(round(20_000_000-PANEL_M*H/2))
if a0<=H: a0=20_000_000
Aref=np.array([a0+i*H for i in range(PANEL_M)],dtype=np.int64); Bref=Aref+H
reference=fit_panel(Aref,Bref)
ref_c=reference['contribution']/reference['contribution'].sum()
ref_ratios=reference['scales']*reference['cell']/reference['median_H']
tr=[]
for panel_id,(A,B,t) in details.items():
    ratios=t['scales']*t['cell']/t['median_H']; c=np.zeros(len(t['scales']))
    for rr,cc in zip(ref_ratios,ref_c):
        if cc<=1e-12: continue
        k=int(np.argmin(np.abs(np.log(ratios/rr)))); c[k]+=cc
    vv=np.divide(c,t['basis'][0],out=np.zeros_like(c),where=t['basis'][0]>0)
    cov=t['basis']@vv; local_rho=cov[1:]/cov[0]
    alpha=math.log(R_REF)/math.log(t['X_over_H']); pred=alpha*local_rho
    abest=float(np.dot(local_rho,t['theory'])/np.dot(local_rho,local_rho)); best=abest*local_rho
    tr.append({'panel_id':panel_id,'target_X':t['median_X'],'target_H':t['median_H'],
               'target_X_over_H':t['X_over_H'],'alpha_theory_log_ratio':alpha,
               'alpha_best_scalar':abest,'alpha_relative_error':(abest-alpha)/alpha,
               'rho1_pred':float(pred[0]),'rho1_target':float(t['theory'][0]),
               'kernel_rmse_theory_alpha':float(np.sqrt(np.mean((pred-t['theory'])**2))),
               'shape_rmse_theory_alpha':float(np.sqrt(np.mean((pred/pred[0]-t['theory']/t['theory'][0])**2))),
               'kernel_rmse_best_alpha':float(np.sqrt(np.mean((best-t['theory'])**2)))})
transfer=pd.DataFrame(tr)
transfer.to_csv(ROOT/'wide_X_universal_spectrum_transfer.csv',index=False,encoding='utf-8-sig')

summary=pd.DataFrame([
    ['X_min',float(panels['median_X'].min())],['X_max',float(panels['median_X'].max())],
    ['orders_of_magnitude_X',float(math.log10(panels['median_X'].max()/panels['median_X'].min()))],
    ['number_of_panels',float(len(panels))],['panels_per_target',3.0],
    ['mean_fraction_L_over_H_lt_4',float(panels['L_over_H_lt_4'].mean())],
    ['sd_fraction_L_over_H_lt_4',float(panels['L_over_H_lt_4'].std(ddof=1))],
    ['mean_fraction_L_over_H_ge_8',float(panels['L_over_H_ge_8'].mean())],
    ['sd_fraction_L_over_H_ge_8',float(panels['L_over_H_ge_8'].std(ddof=1))],
    ['rho1_log_mean',float(amplitude['rho1_times_log_X_over_H'].mean())],
    ['rho1_log_sd',float(amplitude['rho1_times_log_X_over_H'].std(ddof=1))],
    ['minus_log_2',-math.log(2)],
    ['universal_fine_spectrum_transfer_mean_RMSE',float(transfer['kernel_rmse_theory_alpha'].mean())],
    ['universal_fine_spectrum_transfer_min_RMSE',float(transfer['kernel_rmse_theory_alpha'].min())],
],columns=['metric','value'])
summary.to_csv(ROOT/'wide_X_multiscale_scaling_summary.csv',index=False,encoding='utf-8-sig')

print(panels.to_string(index=False)); print(regression.to_string(index=False)); print(summary.to_string(index=False))
