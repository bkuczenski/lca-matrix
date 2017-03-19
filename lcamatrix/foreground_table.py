import re

TAB_LF = '\\\\ \n'


def tex_sanitize(tex):
    tex = re.sub('%', '\\%', tex)
    tex = re.sub('&', '\\&', tex)
    # tex = re.sub('_', '\\\\textunderscore', tex)  # this doesn't work bc filenames have underscores
    return tex


class ForegroundTable(object):
    def __init__(self, fragment, max_cols=8):
        self._f = fragment
        self._pf = fragment.product_flow
        self.max_cols = max_cols

    @property
    def pdim(self):
        return self._f.pdim

    def _table_start(self, aggregate):
        # begin table
        table = '\\subsection{%s}\n%s\n\n{\small from %s}\n' % (tex_sanitize(self._pf.flow['Name']),
                                                                tex_sanitize('; '.join(self._pf.flow['Compartment'])),
                                                                tex_sanitize(str(self._pf.process)))
        fg_cols = min(self.max_cols + 1, self.pdim)
        table += '\n{\\scriptsize\\sffamily\n\\begin{tabularx}{\\textwidth}{|X|%s|' % (
            'c@{~}' * fg_cols)
        if aggregate:
            table += 'c|'
        table += '}\n\\hline\n'
        return table

    @staticmethod
    def _table_end():
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
            if i >= self.max_cols:
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
            if i >= self.max_cols:
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
        xtilde = self._f.x_tilde()
        for row, fg in enumerate(self._f.foreground):
            table += tex_sanitize('(%d) %s [%s]' % (row, fg.flow['Name'], fg.process['SpatialScope']))
            if row >= self.max_cols:
                agg_add = ' & %4.3g' % xtilde[row]
            else:
                agg_add = ' & '
            for i, val in self._fAf.loc[fg].iteritems():
                value = '%4.3g' % val
                if i >= self.max_cols:
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
        xtilde = self._f.x_tilde()
        table = 'Foreground Node Weights $\\tilde{x}$'
        for i, rows in enumerate(xtilde):
            if i >= self.max_cols:
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
                if i >= self.max_cols:
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
            agg = self._f.show_ad_tilde()
            agg_string = '$\\tilde{a_d}$'
        else:
            agg = None
            agg_string = None
        table = self._table_header('Background Dependencies', aggregate=agg_string)
        table += self._do_dep_table(self._f.show_Ad(), agg, max_rows)

        return table

    def _bf_table(self, aggregate=False, max_rows=20):
        if aggregate:
            agg = self._f.show_bf_tilde()
            agg_string = '$\\tilde{b_f}$'
        else:
            agg = None
            agg_string = None
        table = self._table_header('Foreground Emissions and Cutoffs', aggregate=agg_string)
        table += self._do_dep_table(self._f.show_Bf(), agg, max_rows)

        return table

    def foreground_table(self, ad_rows=30, bf_rows=30, max_rows=42, aggregate=False):
        total_rows = len(self._f.foreground)
        table = self._table_start(aggregate)
        if total_rows > max_rows:
            table += self._table_end()
            table += 'Very large foreground ($p=%d$) omitted.\n' % total_rows
            return table
        table += self._af_table(aggregate)
        if aggregate:
            table += self._x_tilde_table()
        total_rows += min(ad_rows, len(self._f.show_Ad()))
        if total_rows > max_rows:
            ad_rows -= (total_rows - max_rows)
            total_rows = max_rows
            bf_rows = 0
        table += self._ad_table(aggregate, ad_rows)
        total_rows += min(bf_rows, len(self._f.show_Bf()))
        if total_rows > max_rows:
            bf_rows -= (total_rows - max_rows)
        table += self._bf_table(aggregate, bf_rows)
        table += self._table_end()
        return table
