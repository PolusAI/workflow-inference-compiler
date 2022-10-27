def main(input_tpr_path, input_trr_path, output_png_path):
    import MDAnalysis as mda
    import nmrformd
    from matplotlib import pyplot as plt

    # The following code comes directly from the nmr4md tutorial at
    # https://github.com/simongravelle/nmrformd/blob/main/docs/source/tutorials/tutorial01.ipynb

    u = mda.Universe(input_tpr_path, input_trr_path)
    nmr_result = nmrformd.NMR(u, ["type H", "type H"])

    fontsize = 20
    font = {'color':  'black', 'weight': 'normal', 'size': fontsize}
    plt.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
        "font.serif": ["Palatino"],
    })

    fig = plt.figure(figsize=(14,7))
    ax1 = plt.subplot(1, 2, 1)
    ax1.loglog(nmr_result.f, 1/nmr_result.R1, 'o', markersize=8) # [:-250]
    ax1.set_xlabel(r"$f$ (MHz)", fontdict=font)
    ax1.set_ylabel(r'$T_1$ (s)', fontdict=font)
    ax1.spines["top"].set_linewidth(2)
    ax1.spines["bottom"].set_linewidth(2)
    ax1.spines["left"].set_linewidth(2)
    ax1.spines["right"].set_linewidth(2)
    ax1.tick_params(axis='x', which='major', pad=10)
    ax1.tick_params(axis='y', which='major', pad=10)
    ax1.minorticks_on()
    ax1.tick_params('both', length=10, width=2, which='major', direction='in')
    ax1.tick_params('both', length=6, width=1.4, which='minor', direction='in')
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    #plt.xlim(5e2, 5e5)
    #plt.ylim(1, 100)
    ax2 = plt.subplot(1, 2, 2)
    ax2.semilogx(nmr_result.t, nmr_result.gij[0], 'o', markersize=8) # [:-250]
    ax2.set_xlabel(r"$t$ (ps)", fontdict=font)
    ax2.set_ylabel(r'$C$', fontdict=font)
    ax2.spines["top"].set_linewidth(2)
    ax2.spines["bottom"].set_linewidth(2)
    ax2.spines["left"].set_linewidth(2)
    ax2.spines["right"].set_linewidth(2)
    ax2.tick_params(axis='x', which='major', pad=10)
    ax2.tick_params(axis='y', which='major', pad=10)
    ax2.minorticks_on()
    ax2.tick_params('both', length=10, width=2, which='major', direction='in')
    ax2.tick_params('both', length=6, width=1.4, which='minor', direction='in')
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    #plt.xlim(5e-1, 5e2)
    #plt.ylim(-0.5e10, 5e10)
    fig.tight_layout()
    plt.savefig(output_png_path)


from workflow_types import *
# NOTE: No other top-level imports supported

inputs = {'input_tpr_path': tprfile,
          'input_trr_path': trrfile,
          'output_png_path': {**string, 'default': 'nmr.png'}}
outputs = {'output_png_path': ('$(inputs.output_png_path)', pngfile)}
