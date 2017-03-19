import pandas as pd
import numpy as np
from scipy.sparse import vstack

import re
# from lcatools.foreground.report import tex_sanitize

TAB_LF = '\\\\ \n'

MAX_COLS = 8  # display node weights in agg column


def tex_sanitize(tex):
    tex = re.sub('%', '\\%', tex)
    tex = re.sub('&', '\\&', tex)
    # tex = re.sub('_', '\\\\textunderscore', tex)  # this doesn't work bc filenames have underscores
    return tex


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

    def _table_start(self, aggregate):
        # begin table
        table = '\\subsection{%s}\n%s\n\n{\small from %s}\n' % (tex_sanitize(self._pf.flow['Name']),
                                                                tex_sanitize('; '.join(self._pf.flow['Compartment'])),
                                                                tex_sanitize(str(self._pf.process)))
        fg_cols = min(MAX_COLS + 1, self.pdim)
        table += '\n{\\scriptsize\\sffamily\n\\begin{tabularx}{\\textwidth}{|X|%s|' % (
            'c@{~}' * fg_cols)
        if aggregate:
            table += 'c|'
        table += '}\n\\hline\n'
        return table

    def _table_end(self):
        return '\\end{tabularx}\n}\n'

    def _table_header(self, title, aggregate=None):
        """
        Add a strut for good measure.
        :param title:
        :param aggregate:
        :return:
        """
        table = title + ' \\rule[-3pt]{0pt}{12pt}'
        for i in range(self.pdim):
            if i >= MAX_COLS:
                table += ' & $\\ldots$'
                break
            table += ' & %d' % i
        if aggregate is not None:
            table += ' & %s' % aggregate
        table += TAB_LF
        table += '\\hline\n'
        return table

    def _table_ellipsis(self, num):
        table = '$\ldots$ (%d rows omitted)' % num
        for i in range(self.pdim):
            table += ' & '
            if i >= MAX_COLS:
                break
        table += TAB_LF
        return table

    def _af_table(self, aggregate):
        # foreground heading
        if aggregate:
            agg_string = ''
        else:
            agg_string = None
        table = self._table_header('(node) Foreground flow', aggregate=agg_string)
        xtilde = self.x_tilde()
        for row, fg in enumerate(self._foreground):
            table += tex_sanitize('(%d) %s [%s]' % (row, fg.flow['Name'], fg.process['SpatialScope']))
            if row >= MAX_COLS:
                agg_add = ' & %4.3g' % xtilde[row]
            else:
                agg_add = ' & '
            for i, val in self.Af.loc[fg].iteritems():
                value = '%4.3g' % val
                if i >= MAX_COLS:
                    table += ' & $\\ldots$'
                    break
                elif i == row:
                    table += ' & \\refbox '
                elif val != 0:
                    table += ' & %s' % value
                else:
                    table += ' & '
            if aggregate:
                table += agg_add
            table += TAB_LF
        table += '\\hline\n'

        return table

    def _x_tilde_table(self):
        xtilde = self.x_tilde()
        table = 'Foreground Node Weights $\\tilde{x}$'
        for i, rows in enumerate(xtilde):
            if i >= MAX_COLS:
                table += ' & $\\ldots$'
                break
            table += ' & %4.3g' % rows
        table += '& '
        table += TAB_LF
        table += '\\hline\n'
        return table

    def _do_dep_table(self, dataframe, agg, max_rows):
        num_rows = 0
        table = ''
        for row, series in dataframe.iterrows():
            num_rows += 1
            if num_rows > max_rows:
                table += self._table_ellipsis(len(dataframe) - max_rows)
                break
            table += '%s' % tex_sanitize(row.table_label())
            for i, val in enumerate(series):
                if i >= MAX_COLS:
                    table += ' & $\\ldots$ '
                    break
                if val == 0:
                    table += ' & '
                else:
                    table += ' & \\dependency'
            if agg is not None:
                table += ' & %5.3g' % agg.loc[row][0]
            table += TAB_LF
        table += '\\hline\n'
        return table

    def _ad_table(self, aggregate=False, max_rows=20):
        if aggregate:
            agg = self.show_ad_tilde()
            agg_string = '$\\tilde{a_d}$'
        else:
            agg = None
            agg_string = None
        table = self._table_header('Background Dependencies', aggregate=agg_string)
        table += self._do_dep_table(self.show_Ad(), agg, max_rows)

        return table

    def _bf_table(self, aggregate=False, max_rows=20):
        if aggregate:
            agg = self.show_bf_tilde()
            agg_string = '$\\tilde{b_f}$'
        else:
            agg = None
            agg_string = None
        table = self._table_header('Foreground Emissions and Cutoffs', aggregate=agg_string)
        table += self._do_dep_table(self.show_Bf(), agg, max_rows)

        return table

    def foreground_table(self, ad_rows=30, bf_rows=30, max_rows=42, aggregate=False):
        total_rows = len(self._foreground)
        table = self._table_start(aggregate)
        if total_rows > max_rows:
            table += self._table_end()
            table += 'Very large foreground ($p=%d$) omitted.\n' % total_rows
            return table
        table += self._af_table(aggregate)
        if aggregate:
            table += self._x_tilde_table()
        total_rows += min(ad_rows, len(self.show_Ad()))
        if total_rows > max_rows:
            ad_rows -= (total_rows - max_rows)
            total_rows = max_rows
            bf_rows = 0
        table += self._ad_table(aggregate, ad_rows)
        total_rows += min(bf_rows, len(self.show_Bf()))
        if total_rows > max_rows:
            bf_rows -= (total_rows - max_rows)
        table += self._bf_table(aggregate, bf_rows)
        table += self._table_end()
        return table
