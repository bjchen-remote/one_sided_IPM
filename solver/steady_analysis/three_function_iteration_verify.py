"""Independent mesh transfer and saved-history audit for the physical runs."""
from pathlib import Path
import json
import numpy as np
from three_function_iteration_run import PhysicalModel, backend, OUT, dump


def transfer(label):
    config = json.loads((OUT/f'{label}_config.json').read_text())
    saved = np.load(OUT/f'{label}_primal_state.npz')
    z, ccoef, rates = saved['z'], saved['ccoef'], saved['rates']
    args = config['arguments']
    rows = []
    for nx, ny in [(65, 33), (129, 65), (257, 129)]:
        p, _, _, qfun, _ = backend.source_problem(nx, ny)
        model = PhysicalModel(p, rates, qfun, args['shape'], args['height'], config['power'],
                              free_rates=not args['fixed_rates'], smooth_c=True,
                              fit_power=args.get('free_power', False),
                              tie_power=config['dynamic_theory_exponent'])
        residual, state = model.evaluate(z, cache=True)
        fixed_C = model.E@ccoef
        fixed_Q = state['Q']-state['C'][:, None]+fixed_C[:, None]
        _, independent = p.evaluate(np.r_[fixed_Q.ravel(), rates])
        raw = independent['raw'].reshape(p.shape)
        hs, aa, mm = model.unpack(z)
        xx = np.linspace(-2., 2., 1001)
        ff = model.mother(xx, hs)
        fp = -ff*np.exp(hs(xx))*(2*xx/(1+xx*xx)+np.log1p(xx*xx)*hs(xx, 1))
        geometry = dict(admissible=bool(model.admissible(z)),
                        maximum_h_coefficient=float(np.max(abs(z[:model.ns]))),
                        logA_range=[float(aa(p.y).min()), float(aa(p.y).max())],
                        center_range=[float(mm(p.y).min()), float(mm(p.y).max())],
                        left_mother_derivative_min=float(fp[xx<0].min()),
                        right_mother_derivative_max=float(fp[xx>0].max()))
        row = dict(nx=nx, ny=ny, with_C_reeliminated=state['metrics'],
                   with_entire_profile_frozen=dict(rms=float(np.sqrt(np.mean(raw*raw))),
                       inf=float(abs(raw).max()), wall_rms=float(np.sqrt(np.mean(raw[0]**2))),
                       discrete_gauges=independent['gauges'].tolist()), geometry=geometry)
        rows.append(row)
        print(label, nx, 'fixed profile rms', row['with_entire_profile_frozen']['rms'],
              'refitted C rms', state['metrics']['rms'], flush=True)
        if nx == args['nx']:
            old = json.loads((OUT/f'{label}_primal.json').read_text())['final']['rms']
            assert abs(old-state['metrics']['rms']) < 1e-9, 'Backend snapshot did not reproduce earlier run'
    return rows


def audit_histories():
    rows = []
    for path in sorted(OUT.glob('*_dual.json')):
        result = json.loads(path.read_text())
        history = result['history']
        row = dict(label=result['label'], status=result['status'], iterations=result['iterations'],
                   initial_rms=history[0]['rms'], final_rms=history[-1]['rms'],
                   minimum_rms=min(t['rms'] for t in history),
                   final_inf=history[-1]['inf'], final_dual_norm=history[-1]['dual_norm'],
                   actual_tau_range=[min(t['tau'] for t in history if 'tau' in t),
                                     max(t['tau'] for t in history if 'tau' in t)],
                   safeguarded_steps=sum(t.get('safeguard_retries', 0)>0 for t in history),
                   final_step_norm=history[-2].get('step_norm'),
                   residual_near_normal_fraction=result['diagnostic']['normal_residual_fraction'])
        rows.append(row)
    return rows


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
    labels = [('dynamic_65_dual', '65x33, tau=1'),
              ('dynamic_65_tau01_dual', '65x33, tau=0.1'),
              ('rich_dynamic_129_dual', '129x65, richer functions'),
              ('rich_freepower_129_dual', '129x65, free width exponent')]
    for label, legend in labels:
        history = json.loads((OUT/f'{label}.json').read_text())['history']
        it = [t['iteration'] for t in history]
        axes[0].semilogy(it, [t['rms'] for t in history], label=legend)
        axes[1].plot(it, [t['dual_norm'] for t in history], label=legend)
    axes[0].axhline(1e-7, color='gray', linestyle='--', label='root max-norm target')
    axes[0].set(title='Full PDE RMS: no convergence to zero', xlabel='Primal-dual iteration', ylabel='PDE RMS')
    axes[1].set(title='Multiplier keeps accumulating residual', xlabel='Primal-dual iteration', ylabel='Multiplier norm')
    axes[0].legend(fontsize=7, loc='lower left')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(OUT/'convergence.png', dpi=180)


def main():
    result = dict(mesh_transfer={label: transfer(label) for label in
                  ['rich_dynamic_129', 'rich_freepower_129']},
                  histories=audit_histories(),
                  limitation='Frozen inherited exterior; no exclusion of exact unified mothers in the half-plane.')
    dump(OUT/'verification.json', result)
    plot()


if __name__ == '__main__':
    main()
