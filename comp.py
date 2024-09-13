import numpy as np
import numba as nb
import networkx as nx
import ctypes as cts
from scipy.sparse import issparse
import os


def get_ss(src_kernels, time_series, time_points, n_boots=500, alpha=0.01, sorting='mergesort', verbose=False):
    """Compute significant sources distribution."""
    n_src, n_ch, n_trials = src_kernels.shape
    n_tp = len(time_points)
    stim = next(idx for idx, tp in enumerate(time_points) if tp > 0)
    j_avg = np.zeros((n_src, n_tp))
    for t in range(n_trials):
        j_avg += np.dot(src_kernels[..., t], time_series[..., t])
    j_avg /= n_trials
    j_baseline_mean = np.expand_dims(np.mean(j_avg[:,:stim], axis=1), axis=1)
    j_avg -= j_baseline_mean
    j_baseline_std = np.expand_dims(np.std(j_avg[:,:stim], axis=1, ddof=1), axis=1)

    boots = np.zeros((n_boots, stim))
    rt = np.random.permutation(n_trials)
    for n in range(n_boots):
        if verbose:
            print(f'Bootstrap n.{n+1}')
        surrogates = np.zeros((n_src, stim))
        for t in rt:
            surrogates += np.dot(src_kernels[:,:,t], time_series[:,np.random.permutation(stim),t])
        boots[n,:] = np.max(np.abs((surrogates/n_trials-j_baseline_mean)/j_baseline_std), axis=0)
    threshold = np.quantile(boots, 1-alpha) * j_baseline_std
    ss = np.abs(j_avg) > threshold
    sorted_src = np.argsort(ss.sum(axis=1), kind=sorting)

    return ss, j_avg, sorted_src, stim


def get_tss(ss, axis=0):
    """Return total significant sources."""
    return ss.sum(axis=axis)


def get_scd(ss, j_avg, axis=0):
    """Return significant current density."""
    return np.sum(ss*np.abs(j_avg), axis=axis)


def get_scs(ss, dist, axis=0):
    """Return significant current scattering averaged over non-zero sources."""
    tss = get_tss(ss, axis=axis)
    #return np.sum(ss*dist, axis=axis) /  (tss + (tss == 0))
    return np.sum(ss*dist, axis=axis)


def get_distances(vert, target, vert_mat=None, measure='euclidean'):
    """Compute distances between a given target and all the source locations."""
    distances = np.zeros((len(vert), 1))
    vert_to_target = np.sqrt(np.sum((vert-target)**2, axis=1))
    if measure == 'euclidean':
        distances = np.expand_dims(vert_to_target, axis=1)
    elif measure == 'geodesic':
        assert vert_mat is not None, "For geodesic distance, vertices connectivity matrix is needed"
        if issparse(vert):
            g = nx.from_scipy_sparse_array(vert_mat)
        else:
            g = nx.from_numpy_array(vert_mat)
        source = np.argmin(vert_to_target) # find the closest source to target coordinates
        path_len = nx.shortest_path_length(g, source=source)
        for i in np.arange(vert.shape[0]):
            if i in path_len.keys():
                distances[i] = path_len[i]
    else:
        raise Exception("Unknown measure - use either 'euclidean' or 'geodesic'")
    return distances


def get_pci(ss_sorted, alg='numba', verbose=False):
    """Return perturbation complexity index."""
    l1,l2 = ss_sorted.shape

    if alg == 'numba':
        comp_series = lz_complexity(ss_sorted, verbose)
    elif alg == 'c':
        dll_name = 'lzc.so'
        if not os.path.isfile(dll_name):
            raise Exception(f'{dll_name} was not found - maybe still needs to be compiled?')
        ss_int = np.ascontiguousarray(ss_sorted, dtype=np.int32)
        dll = cts.CDLL(dll_name)
        pci_func = dll.lz_complexity
        pci_func.argtypes = (
            np.ctypeslib.ndpointer(cts.c_int, ndim=2, shape=ss_int.shape, flags='C_CONTIGUOUS'),
            cts.c_size_t,
            cts.c_size_t)
        pci_func.restype = np.ctypeslib.ndpointer(cts.c_int, shape=(ss_int.shape[1]), flags='C_CONTIGUOUS')
        comp_series = pci_func(ss_int, ss_int.shape[0], ss_int.shape[1])
    else:
        raise Exception("Unknown algorithm - use either 'numba' or 'c'")

    p1 = np.sum(ss_sorted) / (l1 * l2)
    h = -p1 * np.log2(p1) - (1 - p1) * np.log2(1 - p1)
    pci = comp_series * np.log2(l1 * l2) / (l1 * l2 * h)

    return pci


def lz_complexity(ss, verbose=False):
    """Compute Lempel-Ziv complexity."""
    c = 1
    r = 1
    q = 1
    k = 1
    i = 1

    l1, l2 = ss.shape
    clist = []

    while r <= l2:
        if q == r:
            a = i + k - 1
        else:
            a = l1

        found = seek_seq(ss[0:a, q - 1], ss[i:i + k, r - 1])

        if found:
            k += 1
            if i + k > l1:
                clist.append(c)
                r += 1
                if verbose:
                    print(f'r={r}')
                i = 0
                q = r - 1
                k = 1
        else:
            q -= 1
            if q < 1:
                c += 1
                i += k
                if i + 1 > l1:
                    clist.append(c)
                    r += 1
                    if verbose:
                        print(f'r={r}')
                    i = 0
                    q = r - 1
                    k = 1
                else:
                    q = r
                    k = 1

    clist[-1] += 1

    return np.array(clist)


@nb.jit(fastmath=True)
def seek_seq(data, seq):
    """Search for the presence of a given subsequence within the provided data."""
    for i in range(len(data) - len(seq) + 1):
        for j in range(len(seq)):
            if data[i + j] != seq[j]:
                break
        else:
            return True
    return False
