import pandas as pd

from lcamatrix.foreground import ForegroundFragment

class DisplayFragment(ForegroundFragment):
    """
    Uses Pandas for attractive display and functional manipulation of foreground matrices
    """
    @staticmethod
    def _show_nonzero_rows(dataframe):
        return dataframe.loc[~(dataframe == 0).all(axis=1)]

    def show_Af(self):
        return pd.DataFrame(self._af.todense(), index=[k for k in self._foreground])

    def show_Ad(self):
        ad = pd.DataFrame(self.Ad.todense(), index=[k for k in self.bg_flows],
                          columns=[l.process for l in self._foreground])
        return self._show_nonzero_rows(ad)

    def show_Bf(self):
        bf = pd.DataFrame(self.Bf.todense(), index=[k for k in self.emissions],
                          columns=[l.process for l in self._foreground])
        return self._show_nonzero_rows(bf)

    def show_ad_tilde(self, node=0):
        adt = pd.DataFrame(self.ad_tilde(node), index=[k for k in self.bg_flows])
        return self._show_nonzero_rows(adt)

    def show_bf_tilde(self, node=0):
        bft = pd.DataFrame(self.bf_tilde(node), index=self.emissions)
        return self._show_nonzero_rows(bft)

    def show_E(self):
        e = pd.DataFrame(self.E.T.todense(), index=self.emissions, columns=self._qs)
        return self._show_nonzero_rows(e)
