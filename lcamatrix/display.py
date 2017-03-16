import pandas as pd
import numpy as np

import re
# from lcatools.foreground.report import tex_sanitize

TAB_LF = '\\\\ \n'


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
        return len(self._foreground)

    @property
    def ndim(self):
        return self._bg.tstack.ndim

    @property
    def mdim(self):
        return self._bg.mdim

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

    def _table_start(self, aggregate):
        # begin table
        table = tex_sanitize('{\small %s, from %s}\n' % (self._pf.flow, self._pf.process))
        table += '\n{\\scriptsize\\sffamily\n\\begin{tabularx}{\\textwidth}{|X|%s|' % ('c@{~}' * self.pdim)
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
        table += TAB_LF
        return table

    def _af_table(self, aggregate):
        # foreground heading
        if aggregate:
            agg_string = ''
        else:
            agg_string = None
        table = self._table_header('(node) Foreground flow', aggregate=agg_string)
        if aggregate:
            xtilde = self.x_tilde()
        for row, fg in enumerate(self._foreground):
            table += tex_sanitize('(%d) %s [%s]' % (row, fg.flow['Name'], fg.process['SpatialScope']))
            for i, val in self.Af.loc[fg].iteritems():
                if i == row:
                    table += ' & \\refbox '
                elif val != 0:
                    table += ' & %4.3g' % val
                else:
                    table += ' & '
            if aggregate:
                table += '& '
            table += TAB_LF
        table += '\\hline\n'

        return table

    def _x_tilde_table(self):
        xtilde = self.x_tilde()
        table = 'Foreground Node Weights $\\tilde{x}$'
        for rows in xtilde:
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
            for i, val in series.iteritems():
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
        table = self._table_header('Foreground Emissions', aggregate=agg_string)
        table += self._do_dep_table(self.show_Bf(), agg, max_rows)

        return table

    def foreground_table(self, ad_rows=40, bf_rows=20, aggregate=False):
        table = self._table_start(aggregate)
        table += self._af_table(aggregate)
        if aggregate:
            table += self._x_tilde_table()
        table += self._ad_table(aggregate, ad_rows)
        table += self._bf_table(aggregate, bf_rows)
        table += self._table_end()
        return table
