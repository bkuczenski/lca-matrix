import numpy as np
from scipy.sparse import vstack
import uuid

# from lcatools.foreground.report import tex_sanitize
from lcatools.lcia_results import LciaResult, LciaResults


class ForegroundFragment(object):
    """
    A portion of a background archive that can be represented separately from the background.  ForegroundFragments
    are the basis for LCI queries against background archives; are equivalent to performing an LCI query against
    a background archive.

    Requires the background provider to meet the following interface:
     bg.pdim - number of foreground flows
     bg.ndim - number of background flows
     bg.mdim - number of emissions
     bg.is_background(pf)
     bg.foreground(pf) - returns an ordered list of indices (w.r.t. p) downstream of the named product flow (inclusive)
     bg.make_foreground(pf) - returns Af, Ad, Bf (sparse) for product flow

     bg.background_flows() - generates ProductFlows in the background
     bg.emissions - m-list of Emission objects
     bg.construct_sparse(entries, nrows, ncols) - where entries is [[row index, colum index, data]..] - static

    for LCIA:
     bg.compute_bg_lci(ad) - iteratively calculate x, bx for n-dim input vector ad
     bg.compute_lci(pf) - calculate x, bx, bf_tilde for product flow pf
    """
    def __init__(self, bg, flowdb, product_flow):
        """
        instantiates a foreground fragment
        :param bg: a background manager
        :param flowdb: required for compartments and for characterization
        :param product_flow: a ProductFlow known to the background
        """

        self._bg = bg
        self._db = flowdb
        self._pf = product_flow
        self._uuid = uuid.uuid4()
        self._lcia = []  # array of sparse e vectors
        self._qs = []  # array quantities corresponding to rows in lcia

        if bg.is_background(product_flow):
            self._foreground = [product_flow]
        else:
            self._foreground = bg.foreground(product_flow)
        self._af, self._ad, self._bf = bg.make_foreground(product_flow)

        self._is_elem = np.array([self._db.compartments.is_elementary(f.flow) for f in self.emissions])

        self._bx = None  # cached background LCI

        print('Fragment with %d foreground flows' % self.pdim)
        print(' Ad: %dx%d, %d nonzero' % (self.ndim, self.pdim, self._ad.nnz))
        print(' Bf: %dx%d, %d nonzero' % (self.mdim, self.pdim, self._bf.nnz))

    @property
    def uuid(self):
        return str(self._uuid)

    @property
    def foreground(self):
        return self._foreground

    def __len__(self):
        return len(self._foreground)

    @property
    def product_flow(self):
        return self._pf

    @property
    def fg_flows(self):
        for k in self._foreground:
            yield k

    @property
    def bg_flows(self):
        for k in self._bg.background_flows():
            yield k

    @property
    def emissions(self):
        for k in self._bg.emissions:
            yield k

    @property
    def elementary(self):
        for i, k in enumerate(self._bg.emissions):
            if self._is_elem[i]:
                yield k

    @property
    def cutoffs(self):
        for i, k in enumerate(self._bg.emissions):
            if not self._is_elem[i]:
                yield k

    """ Some tension here about where cutoffs belong and how to expose them.
    Background: cutoffs are mathematically identical to emissions, and they cannot be distinguished without access to
    a database of compartments that can tell whether a flow is elementary or not.  That flow database is not provided
    to the BackgroundManager but it is provided here; thus, it is the fragment's job to distinguish cutoffs from
    elementary flows.  (both are called 'emissions', which also includes resource requirements, but could alternately
    be called 'exterior' flows or somesuch).

    Two methods are implemented for finding cutoffs:

    1- the ForegroundFragment provides access to five foreground matrices:
     * Af - foreground matrix
     * Ad - background dependencies
     * Bf - foreground emissions
     * Bf_elementary - Bf sliced to select rows of elementary Emissions
     * Bf_curoff - Bf sliced to select rows of non-elementary Emissions or cutoffs

    It also has five corresponding generators:
     * fg_flows yields ProductFlows in foreground, corresponding to row/column of Af
     * bg_flows yields ProductFlows corresponding to rows of Ad
     * emissions yields Emissions corresponding to rows of Bf
     * elementary yields Emissions whose flow compartments are elementary; also corresponding to rows of Bf_elementary
     * cutoffs yields Emissions whose flow compartments are not elementary; also corresponding to rows of Bf_cutoff

    Alternately,
    2- it has an is_elem property, which returns an np array of boolean values that correspond to entries in emissions.
    This can be used in place of Bf_elementary and Bf_cutoff to slice Bf directly.
    """

    @property
    def is_elem(self):
        return self._is_elem

    @property
    def pdim(self):
        return len(self._foreground)

    @property
    def ndim(self):
        return self._ad.shape[0]

    @property
    def mdim(self):
        return self._bf.shape[0]

    @property
    def tdim(self):
        return len(self._qs)

    @property
    def lcia_methods(self):
        return self._qs

    @property
    def Af(self):
        return self._af

    @property
    def Ad(self):
        return self._ad

    @property
    def Bf(self):
        return self._bf

    @property
    def Bf_cutoff(self):
        return self._bf[~self._is_elem]

    @property
    def Bf_elementary(self):
        return self._bf[self._is_elem]

    def x_tilde(self, node=0):
        if self._foreground is None:
            return np.matrix([[1]])
        return np.linalg.inv(np.eye(self._af.shape[0]) - self._af.todense())[:, node]

    def ad_tilde(self, node=0):
        return self._ad.todense() * self.x_tilde(node)

    def bf_tilde(self, node=0):
        return self._bf.todense() * self.x_tilde(node)

    @property
    def E(self):
        if self.tdim > 0:
            return vstack(self._lcia)
        else:
            return np.matrix([])

    def characterize(self, quantity):
        """
        Generate an e vector for the given quantity using the flow database specified at initialization.
        :param quantity:
        :return:
        """
        if not quantity.is_lcia_method():
            print('Quantity is not an LCIA method.')
            return
        if quantity in self._qs:
            return
        nums = []
        for m, em in enumerate(self._bg.emissions):
            if em.flow.has_characterization(quantity):
                nums.append((0, m, em.flow.cf(quantity)))
            else:
                cf = self._db.lookup_single_cf(em.flow, quantity)
                if cf is not None:
                    em.flow.add_characterization(cf)
                    nums.append((0, m, cf.value))
        e = self._bg.construct_sparse(np.array(nums), 1, self.mdim)
        self._lcia.append(e)
        self._qs.append(quantity)

    def compute_lcia(self, inv):
        if self.tdim == 0:
            return np.array([])
        return self.E * inv

    def lcia_results(self, lci=None, **kwargs):
        """
        :param lci: if supplied, don't generate
        :param kwargs:
        :return: LciaResults
        """
        if lci is None:
            lci = self._bg.lci(self.product_flow, **kwargs)
        l = LciaResults(self.product_flow)
        for q in self.lcia_methods:
            r = LciaResult(q)
            r.add_component(self._uuid, self.product_flow)
            for ex in lci:
                if ex.flow.has_characterization(q):
                    r.add_score(self._uuid, ex, ex.flow.factor(q), 'GLO')  # just continue to punt on location
            l.add(r)
        return l

    def lcia(self):
        return self.fg_lcia() + self.bg_lcia()

    def fg_lcia(self):
        return self.compute_lcia(self.bf_tilde())

    def bg_lcia(self):
        if self._bx is None:
            _, bx = self._bg.compute_bg_lci(self.ad_tilde())
            self._bx = bx
        return self.compute_lcia(self._bx)

    def pf_lcia(self, pf):
        _, bx, bftilde = self._bg.compute_lci(pf)
        return self.compute_lcia(bx + bftilde)
