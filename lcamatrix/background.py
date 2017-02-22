"""
Tarjan's strongly connected components algorithm is recursive.  Python doesn't do well with deep recursion, so
ultimately this code will need to be implemented on a more grown-up language.  For now, however, the recursion
limit test that ships with python reported a segfault at a recursion limit exceeding 19100 -- bigger than ecoinvent!
So for the time being we are safe.
may need to use threading to go higher (see http://stackoverflow.com/questions/2917210/)
Validate recursion depth on a given system using PYTHONROOT/Tools/scripts/find_recursionlimit.py
"""
import os
import sys  # for recursion limit

import scipy as sp
from scipy.sparse import csc_matrix  # , csr_matrix

from lcamatrix.tarjan_stack import TarjanStack
from lcamatrix.product_flow import ProductFlow
from lcamatrix.emission import Emission

from collections import defaultdict, namedtuple
from lcatools.exchanges import comp_dir
from lcatools.entities import LcProcess

MAX_SAFE_RECURSION_LIMIT = 18000  # this should be validated using


MatrixEntry = namedtuple("MatrixEntry", ('parent', 'term', 'value'))  # parent = column; term = row; value > 0 => input


DEFAULT_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class NoMatchingReference(Exception):
    pass


class NoAllocation(Exception):
    pass


def resolve_termination(exchange, terms, strategy):
    """
         'cutoff' - call the flow a cutoff and ignore it
         'mix' - create a new "market" process that mixes the inputs
         'first' - take the first match (alphabetically by process name)
         'last' - take the last match (alphabetically by process name)

    :param exchange:
    :param terms:
    :param strategy:
    :return:
    """
    if len(terms) == 1:
        return terms[0]
    elif len(terms) == 0 or strategy == 'cutoff':
        return None
    elif strategy == 'mix':
        p = LcProcess.new('Market for %s' % exchange.flow['Name'], Comment='Auto-generated')
        p.add_exchange(exchange.flow, comp_dir(exchange.direction), value=float(len(terms)))
        p.add_reference(exchange.flow, comp_dir(exchange.direction))
        for t in terms:
            p.add_exchange(exchange.flow, exchange.direction, value=1.0, termination=t.get_uuid())
        return p
    elif strategy == 'first':
        return [t for t in sorted(terms, key=lambda x:x['Name'])][0]
    elif strategy == 'last':
        return [t for t in sorted(terms, key=lambda x:x['Name'])][-1]
    else:
        raise KeyError('Unknown multi-termination strategy %s' % strategy)


'''
def is_elementary(flow):
    """
    in future, this sholud lookup to a standalone compartment manager
    :param flow:
    :return:
    """
    comp = flow['Compartment'][0]
    if comp == 'air' or comp == 'water' or comp == 'soil' or comp == 'natural resource':
        # Ecoinvent + USLCI
        return True
    elif comp == 'resource':
        # USLCI
        return True
    return False
'''


class BackgroundManager(object):
    """
    Class for managing a collection of linked processes as a coherent technology matrix.
    """
    def __init__(self, archive, data_dir=DEFAULT_DATA_DIR):
        """

        :param data_dir:
        """
        self.archive = archive
        self._data_dir = data_dir  # save A / Bx in sparse form? Bx isn't sparse- A may be quicker than the json
        self._lowlinks = dict()  # dict mapping product_flow key to lowlink -- which is a key into TarjanStack.sccs

        self.tstack = TarjanStack()
        self._a_matrix = None
        self._b_matrix = None

        self._terminations = defaultdict(list)
        self._index_archive()  # dict of reference flows to terminating processes.
        self._rec_limit = len(self.archive.processes())
        if self.required_recursion_limit > MAX_SAFE_RECURSION_LIMIT:
            raise EnvironmentError('This database may require too high a recursion limit-- time to learn lisp.')

        self._interior_incoming = []  # hold interior exchanges before adding them to the component graph

        self._interior = []  # entries in the sparse A matrix-- may not need to save these since they are encoded into A
        self._foreground = []  # entries upstream of the background
        self._product_flows = dict()  # maps product_flow key to index
        self._pf_index = []  # maps index to product_flow

        self._cutoff = []  # entries in the sparse B matrix
        self._emissions = dict()  # maps emission key to index
        self._ef_index = []  # maps index to emission

    @property
    def required_recursion_limit(self):
        return self._rec_limit

    def index(self, product_flow):
        return self._product_flows[product_flow.key]

    def product_flow(self, index):
        return self._pf_index[index]

    def _lowlink(self, product_flow):
        return self._lowlinks[product_flow.key]

    def _add_product_flow(self, pf):
        self._product_flows[pf.key] = pf.index
        self._set_lowlink(pf, pf.index)
        self._pf_index.append(pf)
        self.tstack.add_to_stack(pf)

    def _set_lowlink(self, pf, lowlink):
        """
        Sets lowlink to be the lower of the existing lowlink or the supplied lowlink
        :param pf:
        :param lowlink:
        :return:
        """
        if pf.key in self._lowlinks:
            self._lowlinks[pf.key] = min(self._lowlink(pf), lowlink)
        else:
            self._lowlinks[pf.key] = lowlink

    def _check_product_flow(self, flow, termination):
        """
        returns the product flow if it exists, or None if it doesn't
        :param flow:
        :param termination: the process whose reference flow is flow
        :return:
        """
        if termination is None:
            k = (flow.uuid, None)
        else:
            k = (flow.uuid, termination.uuid)
        if k in self._product_flows:
            return self.product_flow(self._product_flows[k])
        else:
            return None

    def _create_product_flow(self, flow, termination):
        index = len(self._pf_index)
        pf = ProductFlow(index, flow, termination)
        self._add_product_flow(pf)
        return pf

    def _add_emission(self, flow, direction):
        key = (flow.uuid, direction)
        if key in self._emissions:
            return self._ef_index[self._emissions[key]]
        else:
            index = len(self._ef_index)
            ef = Emission(index, flow, direction)
            self._emissions[ef.key] = index
            self._ef_index.append(ef)
            return ef

    def _index_archive(self):
        """
        Creates a dict of reference flows known to the archive.  The dict maps (flow, direction) to a list of
        processes which terminate it.
        :return:
        """
        for p in self.archive.processes():
            for rx in p.references():
                self._terminations[(rx.flow.uuid, comp_dir(rx.direction))].append(p)

    def terminate(self, exch, strategy):
        """
        Find the ProductFlow that terminates a given exchange.  If an exchange has an explicit termination, use it.
        Otherwise, consult a local cache; and ask the archive [slow] if the cache is not populated.
        :param exch:
        :param strategy:
        :return:
        """
        if exch.termination is not None:
            term = self.archive[exch.termination]
        else:
            key = (exch.flow.uuid, exch.direction)
            if key in self._terminations:
                terms = self._terminations[key]
                if len(terms) == 1:
                    term = terms[0]
                else:
                    term = resolve_termination(exch, terms, strategy)
            else:
                self._terminations[key].append(None)
                term = None
        return term

    def _construct_b_matrix(self):
        """
        b matrix only includes emissions from background + downstream processes.
        [foreground processes LCI will have to be computed the foreground way]
        :return:
        """
        if self._b_matrix is not None:
            raise ValueError('B matrix already specified!')
        num_bg = sp.array([[co.term.index, self.tstack.bg_dict(co.parent.index), co.value] for co in self._cutoff
                           if self.tstack.is_background(co.parent.index)])
        mdim = len(self._emissions)
        self._b_matrix = csc_matrix((num_bg[:, 2], (num_bg[:, 0], num_bg[:, 1])), shape=(mdim, self.tstack.ndim))
        self._cutoff = [co for co in self._cutoff if not self.tstack.is_background(co.parent.index)]

    def _construct_a_matrix(self):
        ndim = self.tstack.ndim
        num_bg = sp.array([[self.tstack.bg_dict(i.term.index), self.tstack.bg_dict(i.parent.index), i.value]
                           for i in self._interior
                           if self.tstack.is_background(i.parent.index)])
        self._a_matrix = csc_matrix((num_bg[:, 2], (num_bg[:, 0], num_bg[:, 1])), shape=(ndim, ndim))

    def _update_component_graph(self):
        self.tstack.add_to_graph(self._interior_incoming)  # background should now be up to date
        while len(self._interior_incoming) > 0:
            k = self._interior_incoming.pop()
            if self.tstack.is_background(k.parent.index):
                self._interior.append(k)
            else:
                self._foreground.append(k)

        if self._a_matrix is None and self.tstack.background is not None:
            self._construct_a_matrix()
            self._construct_b_matrix()

    def add_all_ref_products(self, multi_term='first', default_allocation=None):
        for p in self.archive.processes():
            for x in p.reference_entity:
                j = self._check_product_flow(x.flow, p)
                if j is None:
                    self._add_ref_product(x.flow, p, multi_term, default_allocation)
        self._update_component_graph()

    def add_ref_product(self, flow, termination, multi_term='first', default_allocation=None):
        """
        Here we are adding a reference product - column of the A + B matrix.  The termination must be supplied.
        :param flow: a product flow
        :param termination: a process that includes the product flow among its reference exchanges (input OR output)
        :param multi_term: ['first'] specify how to handle ambiguous terminations.  Possible answers are:
         'cutoff' - call the flow a cutoff and ignore it
         'mix' - create a new "market" process that mixes the inputs
         'first' - take the first match (alphabetically by process name)
         'last' - take the last match (alphabetically by process name)
         Not currently implemented.
        :param default_allocation: an LcQuantity to use for allocation if unallocated processes are encountered
        :return:
        """
        j = self._check_product_flow(flow, termination)

        if j is None:
            j = self._add_ref_product(flow, termination, multi_term, default_allocation)
            self._update_component_graph()
        return j

    def _add_ref_product(self, flow, termination, multi_term, default_allocation):
        old_recursion_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(self.required_recursion_limit)

        j = self._create_product_flow(flow, termination)
        self._traverse_term_exchanges(j, multi_term, default_allocation)

        sys.setrecursionlimit(old_recursion_limit)
        return j

    def _traverse_term_exchanges(self, parent, multi_term, default_allocation=None):
        """
        Implements the Tarjan traversal
        :param parent: a ProductFlow
        :param default_allocation:
        :return:
        """
        rx = parent.process.find_reference(parent.flow)
        no_alloc = False
        if not parent.process.is_allocated(rx):
            if default_allocation is not None:
                parent.process.allocate_by_quantity(default_allocation)
            else:
                no_alloc = True

        if no_alloc:
            exchs = []
        else:
            exchs = [x for x in parent.process.exchanges()]

        for exch in exchs:  # unallocated exchanges, including reference exchs
            if exch in parent.process.reference_entity:
                continue
            val = exch[rx]
            if val is None or val == 0:
                # don't add zero entries (or descendants) to sparse matrix
                continue
            term = self.terminate(exch, multi_term)
            if term is None:
                # cutoff
                emission = self._add_emission(exch.flow, exch.direction)  # check, create, and add all at once
                self.add_cutoff(parent, emission, val)
                continue

            # interior flow-- enforce normative direction
            if exch.direction == 'Output':
                val *= -1
            i = self._check_product_flow(exch.flow, term)
            if i is None:
                # not visited -- need to visit
                i = self._create_product_flow(exch.flow, term)
                self._traverse_term_exchanges(i, multi_term, default_allocation)
                # carry back lowlink, if lower
                self._set_lowlink(parent, self._lowlink(i))
            elif self.tstack.check_stack(i):
                # visited and currently on stack - carry back index if lower
                self._set_lowlink(parent, self.index(i))
            else:
                # visited, not on stack- nothing to do
                pass
            self.add_interior(parent, i, val)

        # name an SCC if we've found one
        if self._lowlink(parent) == self.index(parent):
            self.tstack.label_scc(self.index(parent), parent.key)

    def add_cutoff(self, parent, emission, val):
        """
        Create an exchange for a cutoff flow (incl. elementary flows)
        :param parent: product flow- B matrix column
        :param emission: emission - B matrix row
        :param val: raw exchange value
        """
        value = val / parent.inbound_ev
        self._cutoff.append(MatrixEntry(parent, emission, value))

    def add_interior(self, parent, term, val):
        """
        Enforces the convention that interior exchanges are inputs; reference flows are outputs; symmetrically to
        inbound_ev determination in ProductFlow constructore

        :param parent: product flow - A matrix column
        :param term: product flow - A matrix row
        :param val: raw (direction-adjusted) exchange value
        :return:
        """
        value = val / parent.inbound_ev
        self._interior_incoming.append(MatrixEntry(parent, term, value))
