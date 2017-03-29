import xlwt
from math import ceil, log10
from collections import defaultdict

from lcamatrix.product_flow import ProductFlow


class ForegroundPublication(object):
    """
    Create an XLS document reporting the contents of the foreground fragment.  For Kuczenski (2017) JIE

    The main payload is a ForegroundFragment, which must support the following interface:
     fragment.uuid - uuid of ForegroundFragment- to represent its reference flow

     fragment.mdim - size of emissions
     fragment.pdim - size of foreground
     fragment.tdim - number of characterized lcia methods

     fragment.Af, (p x p)
     fragment.Ad, (: x p)
     fragment.Bf (m x p)- sparse Af, Ad, Bf matrices

     fragment.lcia_methods - t-list of quantities
     fragment.foreground - p-list of ProductFlows (Af dims, Ad and Bf columns)
     fragment.bg_flows - ordered ProductFlow (Ad row) generator
     fragment.emissions - ordered Emission (Bf row) generator
     fragment.is_elem - np array of boolean values discriminating between elementary and non-elementary rows of Bf

     fragment.x_tilde(i) - dense column vector of node weights for unit output of ith node (default: canonical i=0)

    for lcia reporting:
     fragment.E - sparse E matrix (must be t x m)
     fragment.fg_lcia() - returns E * bf_tilde
     fragment.bg_lcia() - returns E * bx
     fragment.pf_lcia(product_flow) - returns E * b for specified product flow

    """
    def _create_xls(self):
        self._x = xlwt.Workbook()
        self._wid = defaultdict(int)

    def _save_xls(self, filename):
        for k, val in self._wid.items():
            self._x.get_sheet(k[0]).col(k[1]).width = (val + 2) * 256
        self._x.save(filename)

    def __init__(self, fragment, detail=True, lcia=True, private=None):
        self._detail = detail
        self._lcia = lcia and fragment.tdim > 0
        self._uuid = fragment.uuid

        self._af = fragment.Af.tocoo()
        self._ad = fragment.Ad.tocoo()
        self._bf = fragment.Bf.tocoo()
        self._xtilde = fragment.x_tilde()

        self._x = None  # storage spot for workbook in progress
        self._wid = None

        self._lm = fragment.lcia_methods
        self._ff = fragment.foreground

        self._ff_idx = dict((ff, n) for n, ff in enumerate(self._ff))  # maps entity to Af row/column
        self._ad_idx = dict((ad, n) for n, ad in enumerate(fragment.bg_flows))  # maps entity to Ad row
        self._bf_idx = dict((bf, n) for n, bf in enumerate(fragment.emissions))  # maps entity to Bf row

        self._ff_len = ceil(log10(len(self._ff_idx)))
        self._ad_len = ceil(log10(len(self._ad_idx)))
        self._bf_len = ceil(log10(len(self._bf_idx)))

        if self._lcia:
            self._lm_idx = dict((lm, n) for n, lm in enumerate(self._lm))
            self._lm_len = ceil(log10(len(self._lm_idx)))
        else:
            self._lm_idx = dict()
            self._lm_len = 0

        if private is None:
            self._private = []
        else:
            # list of indices to Ad rows to conceal
            self._private = [self._ad_idx[k] for k in private if isinstance(k, ProductFlow) and fragment.is_bg(k)]
            # (future: also give the option to conceal Af columns)

        ad_tilde = self._ad * self._xtilde
        bf_tilde = self._bf * self._xtilde

        self._ad_seen = [k for i, k in enumerate(fragment.bg_flows) if ad_tilde[i] != 0 and i not in self._private]
        self._bf_seen = [k for i, k in enumerate(fragment.emissions) if bf_tilde[i] != 0 and fragment.is_elem[i]]
        self._co_seen = [k for i, k in enumerate(fragment.emissions) if bf_tilde[i] != 0 and not fragment.is_elem[i]]

        self._scores = dict()
        self._e = None

        if self._lcia:
            # # outmoded: select only characterizations for foreground emissions
            # self._e = fragment.E[:, bf_tilde.nonzero()[0]].tocoo()
            self._e = fragment.E.tocoo()
            self._compute_scores(fragment)

    def _lm_key(self, idx):
        return 'LM%0*d' % (self._lm_len, idx)

    def _ff_key(self, idx):
        return 'FF%0*d' % (self._ff_len, idx)

    def _ad_key(self, idx):
        return 'AD%0*d' % (self._ad_len, idx)

    def _bf_key(self, idx):
        return 'EM%0*d' % (self._bf_len, idx)

    def key(self, entity):
        if entity in self._lm_idx:
            return self._lm_key(self._lm_idx[entity])
        elif entity in self._ff_idx:
            return self._ff_key(self._ff_idx[entity])
        elif entity in self._ad_idx:
            return self._ad_key(self._ad_idx[entity])
        else:
            # falls through to KeyError
            return self._bf_key(self._bf_idx[entity])

    def _compute_scores(self, fragment):
        """
        The scores dict contains the following keys:
        s_tilde: total
        sf_tilde: foreground total
        sx_tilde: bg total
        sx_priv: private bg total
        AD#####: unit impacts for entries in ad_seen

        and the values satisfy:
        s_tilde = sf_tilde + sx_tilde
        sx_tilde = sx_priv + [AD#####' * ad_tilde]  where * indicates dot product

        :param fragment: the ForegroundFragment, which has the LCIA engine
        :return:
        """
        if self._lcia is False:
            return
        self._scores['sf_tilde'] = fragment.fg_lcia()
        self._scores['sx_tilde'] = fragment.bg_lcia().todense()

        sx_priv = None
        ad_tilde = self._ad * self._xtilde
        for i, k in enumerate(fragment.bg_flows):
            if i in self._private:
                _priv = (fragment.pf_lcia(k) * ad_tilde[i]).todense()
                if sx_priv is None:
                    sx_priv = _priv
                else:
                    sx_priv += _priv
                print('x')
            elif k in self._ad_seen:
                self._scores[self.key(k)] = fragment.pf_lcia(k).todense()
                print(self.key(k))
        self._scores['sx_priv'] = sx_priv
        self._scores['s_tilde'] = self._scores['sx_tilde'] + self._scores['sf_tilde']

    def _print_lm(self, lm=None):
        if lm is None:
            return 'Key', 'Origin', 'Identifier', 'ReferenceUnit', 'Descriptor'
        else:
            return self.key(lm), lm.origin, lm.get_external_ref(), lm.unit(), lm['Name']

    def _print_ff(self, ff=None):
        if ff is None:
            return 'Key', 'Origin', 'Identifier', 'ReferenceUnit', 'Descriptor', 'FlowDirection', 'FlowName'
        else:
            if self._ff_idx[ff] == 0:
                _id = self._uuid
            else:
                _id = ''
            return self.key(ff), 'Foreground', _id, ff.flow.unit(), str(ff.process), ff.direction, ff.flow['Name']

    def _print_ad(self, ad=None):
        if ad is None:
            return 'Key', 'Origin', 'Identifier', 'ReferenceUnit', 'Descriptor', 'FlowDirection', 'ReferenceFlow'
        else:
            return (self.key(ad), ad.process.origin, ad.process.get_external_ref(), ad.flow.unit(),
                    str(ad.process), ad.direction, ad.flow.get_external_ref())

    def _print_bf(self, bf=None):
        if bf is None:
            return 'Key', 'Origin', 'Identifier', 'ReferenceUnit', 'Descriptor', 'FlowDirection', 'Compartment'
        else:
            return (self.key(bf), bf.flow.origin, bf.flow.get_external_ref(), bf.flow.unit(), bf.flow['Name'],
                    bf.direction, '; '.join(filter(None, bf.compartment)))

    def publish(self, filename, full=False):
        """
        Write the workbook
        :return:
        """
        self._create_xls()
        self._write_entity_map()
        if self._lcia:
            self._write_lcia()
            if full:
                self._write_matrix('E.T', 'LciaMethod', 'Emission', self._bf_key, self._lm_key,
                                   self._e.T, full=full)
            else:
                self._write_matrix('E', 'LciaMethod', 'Emission', self._lm_key, self._bf_key,
                                   self._e)
            # self._write_matrix('E', 'LciaMethod', 'Emission', self._lm_key, lambda x: self.key(self._bf_char[x]),
            #                   self._e)

        self._write_matrix('Af', 'ForegroundFlow', 'ForegroundNode', self._ff_key, self._ff_key, self._af, full=full)
        self._write_vector('x_tilde', 'ForegroundNode', self._ff_key, self._xtilde)

        if self._detail:
            self._write_matrix('Ad', 'BackgroundDependency', 'ForegroundNode', self._ad_key, self._ff_key, self._ad,
                               full = full)
        self._write_vector('ad_tilde', 'BackgroundDependency', self._ad_key, self._ad * self._xtilde)

        if self._detail:
            self._write_matrix('Bf', 'Emission', 'ForegroundNode', self._bf_key, self._ff_key, self._bf, full=full)
        self._write_vector('bf_tilde', 'Emission', self._bf_key, self._bf * self._xtilde)

        self._save_xls(filename)

    def _update_width(self, col, data):
        self._wid[col] = max(self._wid[col], len(str(data)))

    def _write_row(self, sheet, row, data):
        for i, d in enumerate(data):
            sheet.write(row, i, d)
            self._update_width((sheet.get_name(), i), d)
        row += 1
        return row

    def _write_block(self, sheet, row, ite, printer):
        row = self._write_row(sheet, row, printer())
        for k in ite:
            row = self._write_row(sheet, row, printer(k))
        row += 1
        return row

    def _write_entity_map(self):
        row = 0
        emap = self._x.add_sheet('EntityMap')

        # first lcia methods
        if self._lcia:
            row = self._write_row(emap, row, ('Characterizaton Quantities', ))
            row = self._write_block(emap, row, self._lm, self._print_lm)

        # next, the foreground
        row = self._write_row(emap, row, ('Foreground Flows',))
        row = self._write_block(emap, row, self._ff, self._print_ff)

        # background
        row = self._write_row(emap, row, ('Background Dependencies',))
        row = self._write_block(emap, row, self._ad_seen, self._print_ad)

        # cutoffs-
        row = self._write_row(emap, row, ('Cutoffs',))
        row = self._write_block(emap, row, self._co_seen, self._print_bf)

        row = self._write_row(emap, row, ('Elementary Flows',))
        self._write_block(emap, row, self._bf_seen, self._print_bf)

    def _write_lcia(self):
        scores = self._x.add_sheet('LciaScores')
        row = self._write_row(scores, 0, ('LciaMethod', ) + tuple(self.key(l) for l in self._lm) + ('comment', ))

        def _write_score_line(r, key, comment):
            return self._write_row(scores, r, [key] + [float(k) for k in self._scores[key]] + [comment])

        row = _write_score_line(row, 's_tilde', 'Total LCIA Score')
        row = _write_score_line(row, 'sf_tilde', 'Foreground LCIA')
        row = _write_score_line(row, 'sx_tilde', 'Background LCIA')
        if self._scores['sx_priv'] is not None:
            row = _write_score_line(row, 'sx_priv', 'Private component of background score')
        for ad in self._ad_seen:
            row = _write_score_line(row, self.key(ad), 'Impact score for a unit output of process %s' % self.key(ad))

    def _write_vector(self, sheetname, name, key, array):
        sheet = self._x.add_sheet(sheetname)
        row = self._write_row(sheet, 0, (name, 'Data'))
        for i, k in enumerate(array):
            if k != 0:
                row = self._write_row(sheet, row, (key(i), float(k)))

    def _write_matrix(self, sheetname, rowname, colname, rowkey, colkey, coo, full=False):
        """
        Write out a matrix, either in sparse format (row, col, data) or full format
        :param sheetname:
        :param rowname:
        :param colname:
        :param rowkey:
        :param colkey:
        :param coo:
        :param full:
        :return:
        """
        sheet = self._x.add_sheet(sheetname)
        if full:
            dense = coo.todense()
            row = self._write_row(sheet, 0, [''] + [colkey(i) for i in range(dense.shape[1])])
            for i, data in enumerate(dense):
                if data.any():
                    rowdat = [rowkey(i)]
                    for elem in data.T:
                        if elem == 0.0:
                            rowdat.append('')
                        else:
                            rowdat.append(float(elem))
                    row = self._write_row(sheet, row, rowdat)
        else:
            row = self._write_row(sheet, 0, (rowname, colname, 'Data'))
            rows = coo.row
            cols = coo.col
            data = coo.data
            for i, k in enumerate(data):
                row = self._write_row(sheet, row, (rowkey(rows[i]), colkey(cols[i]), float(k)))


'''
    def publish(self):
        """
        Publish the fragment to the specified XLS file.
        :param detail: [True] publish Ad and Bf as well as ad_tilde and bf_tilde
        :param lcia: [True] If the fragment contains LCIA methods, publish the methods and results
        :param private: [None] include a list of background flows to omit from the publication. If lcia is being
         published, the LCIA scores of the private background processes will be aggregated into an sx_tilde entry.
        :return:
        """
        ff, bd, fe, lm = self._publish_entity_map(private)
        if lcia is True and self._f.tdim > 0:
            self._publish_scores(private)
        self._publish_foreground(detail)
        self._publish_background(detail, private)
        self._publish_emissions(detail)

    def _publish_entity_map(self, private=None):
        """
        Create an entity map for the fragment, concealing background ProductFlows included in private.
        :param private:
        :return: foreground_flows, background_dependencies, foreground_emissions, lcia_methods: lists of ProductFlows,
         Emissions, and quantities as appropriate
        """
        ff = self._f.foreground
        bd = self._f.show_Ad().index
        if private is not None:
            bd = list(set(bd) - set(private))
        bd = sorted(bd, key=lambda x: x.flow['Name'])
        fe = sorted(list(self._f.show_Bf().index), key=lambda x: x.flow)
        lm = self._f.lcia_methods
        self._write_entity_map(ff, bd, fe, lm)
        return ff, bd, fe, lm

    def _write_lcia_methods(self, sheet, row_pointer, ff):
        pass

    def _write_foreground_flows(self, sheet, row_pointer, ff):
        row_pointer = 0
        # header line
        row_pointer = self._write_row(sheet, row_pointer, ('ForegroundFlow', 'Flow', 'Node', 'ReferenceUnit'))
        for h, f in enumerate(ff):
            row_pointer = self._write_row(sheet, row_pointer,
                                          ('FF%d' % h, f.flow['Name'], str(f.process), f.flow.unit()))
        row_pointer += 1
        return row_pointer

    def _write_background_flows(self, sheet, row_pointer, bd):
        row_pointer = self._write_row(sheet, row_pointer, ('BackgroundDependency', 'Source', 'Link', 'Name',
                                                           'ReferenceFlow'))
        for h, b in enumerate(bd):
            row_pointer = self._write_row(sheet, row_pointer,
                                          ('BD%d' % h, b.process.origin, b.process.external_ref, b.process['Name'],
                                           str(b.flow)))
        row_pointer += 1
        return row_pointer
'''