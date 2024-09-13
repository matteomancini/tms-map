import matplotlib.pyplot as plt


def plot_src(ss, stim, bline=100):
    """Plot significant sources map."""
    plt.figure(figsize=(10,10))
    plt.imshow(ss[:, stim-bline:], interpolation='nearest', aspect='auto')


def plot_measures(measure_list, timepoints):
    """Plot sources-based measures."""
    fig, axs = plt.subplots(len(measure_list))
    for m, a in zip(measure_list, axs):
        a.plot(timepoints, m)
    plt.show()


def plot_pci(pci, timepoints):
    """Plot perturbation complexity index over time."""
    fig, ax = plt.subplots()
    ax.plot(timepoints, pci)
    ax.axhline(y = pci[-1], color = 'k', linestyle = '--')
    font_size = 14
    ax.annotate(f'PCI = {pci[-1]:.2f}', xy=(0,pci[-1]),
                xytext=(0,pci[-1]-font_size*2e-3),fontsize=font_size)
    ax.set(xlabel='Time [ms]', ylabel='PCI', title='PCI')
    plt.show()
