"""Render the new data competition and distinguish two residual diagnostics."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from fresh_positive_search import OUT, sources, wall_model


def main():
    result=json.loads((OUT/'results.json').read_text())
    native=next(sources()); x=np.linspace(-2.5,2.5,601)
    q=RectBivariateSpline(native['y'],native['x'],native['q'])(0,x)[0]
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax[0,0].plot(x,q,color='black',lw=2.5,label='Native data, tau=20')
    labels={'student':'Symmetric power','log_skew':'Log skew','algebraic_warp':'Algebraic warp','atan_skew':'Arctan skew'}
    colors={'student':'#8c8c8c','log_skew':'#b46b28','algebraic_warp':'#156caa','atan_skew':'#4a9476'}
    rows=result['sources']['tau20']['wall_models']
    for row in rows:
        name=row['family']; qh,_=wall_model(x,row['parameters'],name)
        ax[0,0].plot(x,qh,label=labels[name],color=colors[name],lw=1.5)
        ax[0,1].plot(x,qh-q,color=colors[name],label=labels[name])
    for a in ax[0]:
        a.axvspan(-1.25,1.25,color='#dedede',alpha=.4,zorder=-2)
        a.set_xlabel('Normalized physical x')
    ax[0,0].set(ylabel='q = Q_x',title='Two shape parameters suffice for a close fit')
    ax[0,0].legend(fontsize=8)
    ax[0,1].axhline(0,color='black',lw=.5)
    ax[0,1].set(ylabel='Model minus measured q',title='Grey region: fit; outer wings: spatial validation')
    for name,color in [('tau20','#156caa'),('L4_257x129','#b46b28')]:
        rr=result['sources'][name]['native_selected_shape_audit']['heights']
        y=[r['y'] for r in rr]
        ax[1,0].semilogy(y,[r['residual_term_relative_l2'] for r in rr],':',color=color,label=name+': free projected derivative')
        ax[1,0].semilogy(y,[r['actual_geometry_transport_term_relative_l2'] for r in rr],'-o',ms=3,color=color,label=name+': derivative of fitted geometry')
    ax[1,0].set(xlabel='Height y',ylabel='Transport residual / term scale',title='Small projected residual is not closure')
    ax[1,0].legend(fontsize=8)
    t=result['temporal_validation']['records']
    ax[1,1].plot([r['tau'] for r in t],[r['refit_parameters'][0] for r in t],'-o',label='p: fitted exponent',color='#156caa')
    ax[1,1].plot([r['tau'] for r in t],[r['refit_parameters'][1] for r in t],'-s',label='k: fitted skew',color='#4a9476')
    ax[1,1].set(xlabel='Native trajectory time tau',ylabel='Refitted shape parameter',title='Finite-time variation remains measurable')
    ax[1,1].legend(fontsize=8)
    for a in ax.flat:
        a.grid(alpha=.15)
        a.spines[['top','right']].set_visible(False)
    fig.suptitle('Fresh positive-skeleton search — numerical fit is not an exact identity',fontsize=14)
    fig.savefig(OUT/'comparison.png',dpi=180)
    fig.savefig(OUT/'comparison.pdf')
    plt.close(fig)


if __name__=='__main__':main()
