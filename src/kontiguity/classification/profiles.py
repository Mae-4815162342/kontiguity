"""Contact profiles of candidate sequences against the host chromosomes, read from a (m)cool file.

For a candidate c and every host bin i (host chromosomes concatenated):
    r_i  = raw trans contacts between c and bin i
    e_i  = R * w_i,  R = sum_i r_i,  w_i = coverage_i / sum(coverage)   (coverage-based expectation, SCN-like)
Host bins whose coverage is < 10 % of the median are masked.
"""
import h5py
import numpy as np
import scipy.sparse as sp


def _group(f, resolution):
    if "resolutions" in f:
        res = sorted(int(r) for r in f["resolutions"])
        if resolution not in res:
            raise ValueError(f"resolution {resolution} not in mcool (available: {res})")
        return f["resolutions"][str(resolution)]
    return f


class Profiles:
    def __init__(self, path, host_names, resolution=5000):
        if "::" in path:
            path, grp = path.split("::")
            resolution = int(grp.rstrip("/").split("/")[-1]) if "resolutions" in grp else resolution
        with h5py.File(path, "r") as f:
            g = _group(f, resolution)
            self.names = np.array([x.decode() if isinstance(x, bytes) else x for x in g["chroms/name"][:]])
            self.lengths = g["chroms/length"][:].astype(np.int64)
            off = g["indexes/chrom_offset"][:].astype(np.int64)
            b1 = g["pixels/bin1_id"][:].astype(np.int64)
            b2 = g["pixels/bin2_id"][:].astype(np.int64)
            c = g["pixels/count"][:].astype(np.float64)
        nb = off[-1]
        ns = len(self.names)
        cob = np.repeat(np.arange(ns), np.diff(off))
        nd = b1 != b2
        cov = np.bincount(b1[nd], c[nd], nb) + np.bincount(b2[nd], c[nd], nb)

        is_host = np.isin(self.names, list(host_names))
        if not is_host.any():
            raise ValueError("no host chromosome found in the contact map")
        host_bins = np.where(is_host[cob])[0]
        hidx = -np.ones(nb, np.int64); hidx[host_bins] = np.arange(len(host_bins))
        hcov = cov[host_bins]
        med = np.median(hcov[hcov > 0]) if (hcov > 0).any() else 0
        self.hmask = hcov > 0.1 * med
        w = np.where(self.hmask, hcov, 0.0)
        self.w = w / w.sum()
        self.chrom = cob[host_bins]
        self.bin_start = (host_bins - off[cob[host_bins]]) * resolution
        self.host = is_host

        cand = np.where(~is_host)[0]
        self.cand_names = self.names[cand]
        self.cand_lengths = self.lengths[cand]
        cidx = -np.ones(ns, np.int64); cidx[cand] = np.arange(len(cand))
        cb = cidx[cob]
        m1 = (cb[b1] >= 0) & (hidx[b2] >= 0)
        m2 = (cb[b2] >= 0) & (hidx[b1] >= 0)
        self.CH = sp.csr_matrix((np.r_[c[m1], c[m2]], (np.r_[cb[b1[m1]], cb[b2[m2]]], np.r_[hidx[b2[m1]], hidx[b1[m2]]])),
                                shape=(len(cand), len(host_bins)))
        self.CH.sum_duplicates()
        self.trans = np.asarray(self.CH.sum(1)).ravel()
        ms = (cb[b1] >= 0) & (cb[b1] == cb[b2])
        self.self_contacts = np.bincount(cb[b1[ms]], c[ms], len(cand))
        # sequence x sequence contact matrix (for the Hi-C communities)
        S = sp.coo_matrix((c, (cob[b1], cob[b2])), shape=(ns, ns)).tocsr()
        self.seq_contacts = (S + S.T - sp.diags(S.diagonal())).tocsr()
        self.index = {n: i for i, n in enumerate(self.cand_names)}

    def counts(self, names):
        idx = [self.index[n] for n in names]
        return np.asarray(self.CH[idx].sum(0)).ravel()
