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
from scipy.sparse import csc_matrix, csr_matrix

from collections import defaultdict, namedtuple
from lcatools.exchanges import comp_dir, NoAllocation
from lcatools.entities import LcProcess

MAX_SAFE_RECURSION_LIMIT = 18000  # this should be validated using


MatrixEntry = namedtuple("MatrixEntry", ('parent', 'term', 'value'))  # parent = column; term = row; value > 0 => input


DEFAULT_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class NoMatchingReference(Exception):
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


class TarjanStack(object):
    """
    Stores the current stack and provides a record of named SCCs
    """
    def __init__(self):
        self._stack = []
        self._stack_hash = set()
        self._sccs = defaultdict(set)  # dict mapping lowest index (lowlink = SCC ID) to its scc peers
        self._scc_of = dict()  # dict mapping product flow to SCC ID (reverse mapping of _sccs)

        self._component_cols_by_row = defaultdict(set)  # dict of scc_id to upstream dependents (columns with nz row)
        self._component_rows_by_col = defaultdict(set)  # rows with nz column

        self._background = None  # single scc_id representing largest scc
        self._downstream = set()  # sccs on which background depends

    def check_stack(self, product_flow):
        """
        :param product_flow:
        :return:
        """
        return product_flow in self._stack_hash

    def add_to_stack(self, product_flow):
        if not isinstance(product_flow, ProductFlow):
            raise TypeError('TarjanStack should consist only of ProductFlows')
        if self.check_stack(product_flow):
            raise ValueError('ProductFlow already in stack')
        self._stack.append(product_flow)
        self._stack_hash.add(product_flow)

    def label_scc(self, index, key):
        """

        :param index: the index of the lowest link in the SCC-- becomes scc ID
        :param key: the identifier for the lowest link in the SCC (necessary to ID the link)
        :return:
        """
        while 1:
            node = self._stack.pop()
            self._stack_hash.remove(node)
            self._sccs[index].add(node)
            self._scc_of[node] = index
            if node.key == key:
                break

    def _set_background(self):
        ml = 0
        ind = None
        for i, t in self._sccs.items():
            if len(t) > ml:
                ml = len(t)
                ind = i
        self._background = ind

    def _set_downstream(self, upstream=None):
        """
        recursive function to tag all nodes downstream of the named node.
        :param upstream: [None] if none, use background
        :return:
        """
        if upstream is None:
            upstream = self._background

        for dep in self._component_rows_by_col[upstream]:
            if dep != upstream:  # skip self-dependencies
                self._downstream.add(dep)
                self._set_downstream(dep)

    def add_to_graph(self, interiors):
        """
        take a list of interior exchanges (parent, term, exch) and add them to the component graph
        :return:
        """
        for i in interiors:
            row = self.scc_id(i.term)
            col = self.scc_id(i.parent)
            self._component_cols_by_row[row].add(col)
            self._component_rows_by_col[col].add(row)
        self._set_background()
        self._set_downstream()

    @property
    def background(self):
        return self._background

    def background_flows(self):
        """
        Generator. Yields product flows in the db background or downstream.
        :return:
        """
        for i in self.scc(self._background):
            yield i
        for i in self._downstream:
            for j in self.scc(i):
                yield j

    def scc_id(self, pf):
        return self._scc_of[pf]

    def scc(self, index):
        return self._sccs[index]

    def scc_peers(self, pf):
        """
        Returns nodes in the same SCC as product flow
        :param pf:
        :return:
        """
        for i in self._sccs[self._scc_of[pf]]:
            yield i


class ProductFlow(object):
    """
    Class for storing foreground-relevant information about a single matched row-and-column in the interior matrix.

    """
    def __init__(self, index, flow, process):
        """
        Initialize a row+column in the technology matrix.  Each row corresponds to a reference exchange in the database,
        and thus represents a particular process generating / consuming a particular flow.  A ProductFlow entry is
        akin to a fragment termination.

        :param flow: the LcFlow entity that represents the commodity (term_flow in the fragment sense)
        :param process: the termination of the parent node's exchange (term_node). None is equivalent to a
        cutoff flow or elementary flow (distinction is left to a compartment manager).  If non-null, the process must
        possess a reference exchange with the same flow or the graph traversal may be curtailed.
        """
        self._index = index
        self._flow = flow
        self._process = process

        self._hash = (flow.get_uuid(), None)
        self._inbound_ev = 1.0

        if process is None:
            raise TypeError('No termination? should be a cutoff.')

        if len([x for x in process.reference_entity if x.flow == flow]) == 0:
            # still a cutoff- raise a flag but not an error
            print('NoMatchingReference: Flow: %s, Termination: %s' % (flow.uuid, process.uuid))
        else:
            self._hash = (flow.uuid, process.uuid)
            ref_exch = process.reference(flow)
            self._inbound_ev = ref_exch.value
            if self._inbound_ev is None:
                print('None inbound ev! using 1.0. f:%s t:%s' % (flow, process))
                self._inbound_ev = 1.0
            if ref_exch.direction == 'Input':
                self._inbound_ev *= -1

    def __eq__(self, other):
        """
        shortcut-- allow comparisons without dummy creation
        :param other:
        :return:
        """
        return hash(self) == hash(other)
        # if not isinstance(other, ProductFlow):
        #    return False
        # return self.flow == other.flow and self.process == other.process

    def __hash__(self):
        return hash(self._hash)

    @property
    def index(self):
        return self._index


    @property
    def key(self):
        return self._hash

    @property
    def flow(self):
        return self._flow

    @property
    def process(self):
        return self._process

    @property
    def inbound_ev(self):
        return self._inbound_ev


class Emission(object):
    """
    Class for storing exchange information about a single row in the exterior (cutoff) matrix.

    """
    def __init__(self, index, flow, direction):
        """
        Initialize a row+column in the technology matrix.  Each row corresponds to a reference exchange in the database,
        and thus represents a particular process generating / consuming a particular flow.  A ProductFlow entry is
        akin to a fragment termination.

        :param flow: the LcFlow entity that represents the commodity (term_flow in the fragment sense)
        :param direction: the direction of the exchange
        """
        self._index = index
        self._flow = flow
        self._direction = direction

        self._hash = (flow.get_uuid(), direction)

    def __eq__(self, other):
        """
        shortcut-- allow comparisons without dummy creation
        :param other:
        :return:
        """
        return hash(self) == hash(other)
        # if not isinstance(other, ProductFlow):
        #    return False
        # return self.flow == other.flow and self.process == other.process

    def __hash__(self):
        return hash(self._hash)

    @property
    def index(self):
        return self._index

    @property
    def key(self):
        return self._hash

    @property
    def flow(self):
        return self._flow

    @property
    def direction(self):
        return self._direction


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
        self._lowlinks = dict()  # dict mapping product_flow to lowlink -- which is a key into TarjanStack.sccs

        self.tstack = TarjanStack()
        self._a_matrix = None
        self._b_matrix = None

        self._terminations = self._index_archive()  # dict of reference flows to terminating processes.

        self._interior_incoming = []  # hold interior exchanges before adding them to the component graph

        self._interior = []  # entries in the sparse A matrix
        self._cutoff = []  # entries in the sparse B matrix
        self._product_flows = dict()  # maps product_flow to index
        self._pf_index = []  # maps index to product_flow

        self._emissions = dict()
        self._ef_index = []

    def index(self, product_flow):
        return self._product_flows[product_flow]

    def product_flow(self, index):
        return self._pf_index[index]

    def _lowlink(self, product_flow):
        return self._lowlinks[product_flow]

    def _add_product_flow(self, pf):
        self._product_flows[pf] = pf.index
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
        if pf in self._lowlinks:
            self._lowlinks[pf] = min(self._lowlink(pf), lowlink)
        else:
            self._lowlinks[pf] = lowlink

    def _check_product_flow(self, flow, termination):
        """
        returns the product flow if it exists, or None if it doesn't
        :param flow:
        :param termination:
        :return:
        """
        if termination is None:
            k = (flow.uuid, None)
        else:
            k = (flow.uuid, termination.uuid)
        if k in self._product_flows:
            return self.product_flow(self.index(k))
        else:
            return None

    def _create_product_flow(self, flow, termination):
        index = len(self._pf_index)
        pf = ProductFlow(index, flow, termination)
        self._add_product_flow(pf)
        return pf

    def _add_emission(self, flow, direction):
        key = (flow, direction)
        if key in self._emissions:
            return self._ef_index[self._emissions[key]]
        else:
            index = len(self._emissions)
            ef = Emission(index, flow, direction)
            self._emissions[ef] = index
            self._ef_index.append(ef)
            return ef

    def _index_archive(self):
        """
        Creates a dict of reference flows known to the archive.  The dict maps (flow, direction) to a list of
        processes which terminate it.
        :return:
        """
        ref_dict = defaultdict(list)
        for p in self.archive.processes():
            for rx in p.references():
                ref_dict[(rx.flow, comp_dir(rx.direction))].append(p)
        return ref_dict

    def terminate(self, exch, strategy):
        """
        Find the ProductFlow that terminates a given exchange.  If an exchange has an explicit termination, use it.
        Otherwise, consult a local cache; and ask the archive [slow] if the cache is not populated.
        :param exch:
        :param strategy:
        :return:
        """
        key = (exch.flow, exch.direction)
        if exch.termination is not None:
            term = self.archive[exch.termination]
        else:
            if key in self._terminations:
                terms = self._terminations[key]
                term = resolve_termination(exch, terms, strategy)
            else:
                term = None
        return term

    def _construct_b_matrix(self, bg_dict):
        """
        bg_dict maps pf index to column index.
        :param bg_dict:
        :return:
        """
        if self._b_matrix is not None:
            raise ValueError('B matrix already specified!')
        num_bg = sp.array([[co.term.index, bg_dict[co.parent.index], co.value] for co in self._cutoff
                           if co.parent.index in bg_dict])
        self._b_matrix = csc_matrix((num_bg[:, 2], (num_bg[:, 0], num_bg[:, 1])))
        self._cutoff = [co for co in self._cutoff if co.parent.index not in bg_dict]

    def _construct_a_matrix(self):
        bg = [b.index for b in self.tstack.background_flows()]  # list of indices to background nodes
        bg_dict = dict((ind, n) for n, ind in enumerate(bg))  # mapping of pf index to a-matrix index / b-col index
        num_bg = sp.array([[bg_dict[i.term.index], bg_dict[i.parent.index], i.value] for i in self._interior
                           if i.parent.index in bg_dict])
        self._a_matrix = csc_matrix((num_bg[:, 2], (num_bg[:, 0], num_bg[:, 1])))
        #self._construct_b_matrix(bg_dict)

    def _update_component_graph(self):
        self.tstack.add_to_graph(self._interior_incoming)  # background should now be up to date
        self._interior.extend(self._interior_incoming)
        self._interior_incoming = []

        if self._a_matrix is None and self.tstack.background is not None:
            self._construct_a_matrix()

    def add_ref_product(self, flow, termination, multi_term='first', default_allocation=None):
        """
        Here we are adding a reference product - column of the A + B matrix.  The termination must be supplied.
        :param flow:
        :param termination:
        :param multi_term: ['first'] specify how to handle ambiguous terminations.  Possible answers are:
         'cutoff' - call the flow a cutoff and ignore it
         'mix' - create a new "market" process that mixes the inputs
         'first' - take the first match (alphabetically by process name)
         'last' - take the last match (alphabetically by process name)
         Not currently implemented.
        :param default_allocation: an LcQuantity to use for allocation if unallocated processes are encountered
        :return:
        """
        old_recursion_limit = sys.getrecursionlimit()
        required_recursion_limit = len(self.archive.processes())
        if required_recursion_limit > MAX_SAFE_RECURSION_LIMIT:
            raise EnvironmentError('This database may require too high a recursion limit-- time to learn lisp.')
        sys.setrecursionlimit(required_recursion_limit)
        j = self._check_product_flow(flow, termination)
        if j is None:
            j = self._create_product_flow(flow, termination)
            self._traverse_term_exchanges(j, multi_term, default_allocation)

        sys.setrecursionlimit(old_recursion_limit)
        self._update_component_graph()
        return j

    def _traverse_term_exchanges(self, parent, multi_term, default_allocation=None):
        """
        Implements the Tarjan traversal
        :param parent: a ProductFlow
        :param default_allocation:
        :return:
        """
        try:
            exchs = [x for x in parent.process.exchanges(parent.flow)]
        except NoAllocation:
            if default_allocation is not None:
                parent.process.allocate_by_quantity(default_allocation)
                exchs = [x for x in parent.process.exchanges(parent.flow)]
            else:
                raise
        for exch in exchs:  # allocated exchanges, excluding reference exchs
            term = self.terminate(exch, multi_term)
            if term is None:
                self.add_cutoff(parent, exch)
                continue
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
            self.add_interior(parent, exch, i)

        # name an SCC if we've found one
        if self._lowlink(parent) == self.index(parent):
            self.tstack.label_scc(self.index(parent), parent.key)

    def add_cutoff(self, parent, exchange):
        """
        Create an exchange for a cutoff flow (incl. elementary flows)
        :param parent:
        :param exchange: an LcExchange belonging to the parent node
        """
        if exchange.value is None or exchange.value == 0:
            # don't add zero flows to a sparse matrix
            return
        emission = self._add_emission(exchange.flow, exchange.direction)
        value = exchange.value / parent.inbound_ev
        self._cutoff.append(MatrixEntry(parent, emission, value))

    def add_interior(self, parent, exchange, term):
        """
        Enforces the convention that interior exchanges are inputs; reference flows are outputs; symmetrically to
        inbound_ev determination in ProductFlow constructore

        :param parent:
        :param exchange:
        :param term:
        :return:
        """
        if exchange.value is None or exchange.value == 0:
            return
        value = exchange.value / parent.inbound_ev
        if exchange.direction == 'Output':
            value *= -1
        self._interior_incoming.append(MatrixEntry(parent, term, value))
