from ..imports import *
from .CommonPanels import *

def physical_summary(pops):
    fi = plt.figure(figsize=(10, 5), dpi=200)
    gs = plt.matplotlib.gridspec.GridSpec(2, 3,
                                          width_ratios=[1, 4, 1],
                                          height_ratios=[1, 3],
                                          wspace=0.04,
                                          hspace=0.04,
                                          bottom=0.21,
                                          top=0.98,
                                          right=0.98,
                                          left=0.08)#, sharey='row', sharex='col')


    fr = FluxRadius()
    fr.build(pops=pops,
             ax=plt.subplot(gs[1,1]))
    fr.plot_hz()
    fr.ticks_simplify_exponents('y')
    fr.add_teq_axis()

    fr.remove_ylabel()

    mr = MassRadius()
    mr.build(pops=pops,
             ax=plt.subplot(gs[1,0], sharey=fr.ax))
    mr.plot_both_seager(zorder=1e9)
    mr.ticks_enforce_multiple_oom('x')
    mr.ticks_simplify_exponents('xy')


    er = FluxEscape()
    er.build(pops=pops,
             ax=plt.subplot(gs[0,1], sharex=fr.ax))
    er.remove_xlabel()

    er.plot_constant_lambda()


    sr = StellarRadius()
    sr.build(pops=pops,
             ax=plt.subplot(gs[1,2], sharey=fr.ax))
    sr.remove_ylabel()
    mr.ticks_simplify_exponents('y')


    label = '-'.join(list(pops.keys()))
    plt.savefig(f'{label}.pdf')
