import pandas as pd
import numpy as np
from scipy.sparse import vstack

# from lcatools.foreground.report import tex_sanitize


class ForegroundFragment(object):
    """
    Class for easy display of foreground data
    """
    def __init__(self, bg, product_flow):
        """
        instantiates a foreground fragment
        :param bg: a background manager
        :param product_flow: a ProductFlow known to the background
        """
        self._bg = bg
        self._pf = product_flow
        self._lcia = []  # array of sparse e vectors
        self._qs = []  # array quantities corresponding to rows in lcia

        if bg.tstack.is_background(product_flow):
            self._foreground = None
            self._af = np.matrix([[0]])
            self._ad, self._bf = bg.make_background(product_flow)
        else:
            self._foreground = bg.tstack.foreground(product_flow)
            self._af, self._ad, self._bf = bg.make_foreground(self._foreground)
        print('Fragment with %d foreground flows' % self.pdim)
        print(' Ad: %dx%d, %d nonzero' % (self.ndim, self.pdim, self._ad.nnz))
        print(' Bf: %dx%d, %d nonzero' % (self.mdim, self.pdim, self._bf.nnz))

    @property
    def foreground(self):
        return self._foreground

    @property
    def product_flow(self):
        return self._pf

    def __len__(self):
        return len(self._foreground)

    @property
    def pdim(self):
        if self._foreground is None:
            return 0
        return len(self._foreground)

    @property
    def ndim(self):
        return self._bg.tstack.ndim

    @property
    def mdim(self):
        return self._bg.mdim

    @property
    def tdim(self):
        return len(self._qs)

    @property
    def lcia_methods(self):
        return self._qs

    @property
    def Af(self):
        return pd.DataFrame(self._af.todense(), index=[k for k in self._foreground])

    @property
    def Ad(self):
        return self._ad

    @property
    def Bf(self):
        return self._bf

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

    def characterize(self, flowdb, quantity):
        """
        Generate an e vector for the given quantity using the specified flow database.
        :param flowdb:
        :param quantity:
        :return:
        """
        if not quantity.is_lcia_method():
            print('Quantity is not an LCIA method.')
            return
        nums = []
        for m, em in enumerate(self._bg.emissions):
            val = flowdb.lookup_cf_from_flowable(em.flow, em.compartment, quantity)
            if val is not None:
                nums.append((0, m, val.value))
        e = self._bg.construct_sparse(np.array(nums), 1, self.mdim)
        self._lcia.append(e)
        self._qs.append(quantity)

    def lcia(self, inv):
        if self.tdim == 0:
            return np.array([])
        return self.E * inv

    def fg_lcia(self):
        return self.lcia(self.bf_tilde())

    def bg_lcia(self):
        _, bx = self._bg.compute_bg_lci(self.ad_tilde())
        return self.lcia(bx)

    @staticmethod
    def _show_nonzero_rows(dataframe):
        return dataframe.loc[~(dataframe == 0).all(axis=1)]

    def show_Ad(self):
        ad = pd.DataFrame(self.Ad.todense(), index=[k for k in self._bg.tstack.background_flows()],
                          columns=[l.process for l in self._foreground])
        return self._show_nonzero_rows(ad)

    def show_Bf(self):
        bf = pd.DataFrame(self.Bf.todense(), index=self._bg.emissions, columns=[l.process for l in self._foreground])
        return self._show_nonzero_rows(bf)

    def show_ad_tilde(self, node=0):
        adt = pd.DataFrame(self.ad_tilde(node),
                           index=[k for k in self._bg.tstack.background_flows()])
        return self._show_nonzero_rows(adt)

    def show_bf_tilde(self, node=0):
        bft = pd.DataFrame(self.bf_tilde(node), index=self._bg.emissions)
        return self._show_nonzero_rows(bft)

    def show_E(self):
        e = pd.DataFrame(self.E.T.todense(), index=self._bg.emissions, columns=self._qs)
        return self._show_nonzero_rows(e)
