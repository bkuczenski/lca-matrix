import os
from collections import defaultdict, namedtuple
from lcatools.exchanges import comp_dir, NoAllocation
from lcatools.entities import LcProcess

MatrixEntry = namedtuple("MatrixEntry", ('parent', 'term', 'exchange'))

DEFAULT_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class NoMatchingReference(Exception):
    pass


def ambiguous_termination(exchange, terms, strategy):
    """
         'cutoff' - call the flow a cutoff and ignore it
         'mix' - create a new "market" process that mixes the inputs
         'first' - take the first match (alphabetically by process name)
         'last' - take the last match (alphabetically by process name)

    :param exchange:
    :param terms: a list of terminations for the exchange
    :param strategy:
    :return:
    """
    if len(terms) == 1:
        return terms[0]
    elif len(terms) == 0 or strategy == 'cutoff':
        return None
    elif strategy == 'mix':
        p = LcProcess.new('Market for %s' % exchange.flow['Name'], Comment='Auto-generated')
        p.add_reference(exchange.flow, comp_dir(exchange.direction))
        p.add_exchange(exchange.flow, comp_dir(exchange.direction), value=float(len(terms)))
        for t in terms:
            p.add_exchange(exchange.flow, exchange.direction, value=1.0, termination=t.get_uuid())
        return p
    elif strategy == 'first':
        return [t for t in sorted(terms, key=lambda x:x['Name'])][0]
    elif strategy == 'last':
        return [t for t in sorted(terms, key=lambda x:x['Name'])][-1]
    else:
        raise KeyError('Unknown multi-termination strategy %s' % strategy)


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


class TarjanStack(object):
    """
    Stores the current stack and provides a record of named SCCs
    """
    def __init__(self):
        self._stack = []
        self._sccs = defaultdict(set)

    def check_stack(self, product_flow):
        """
        :param product_flow:
        :return:
        """
        return product_flow in self._stack

    def add_to_stack(self, product_flow):
        if not isinstance(product_flow, ProductFlow):
            raise TypeError('TarjanStack should consist only of ProductFlows')
        self._stack.append(product_flow)

    def label_scc(self, index):
        """

        :param index:
        :return:
        """
        while 1:
            node = self._stack.pop()
            self._sccs[index].add(node)
            if node.key == index:
                break

    @property
    def sccs(self):
        return self._sccs.keys()

    def scc(self, key):
        for i in self._sccs[key]:
            yield i


class ProductFlow(object):
    """
    Class for storing foreground-relevant information about a single matched row-and-column in the interior matrix.

    """
    def __init__(self, flow, process=None):
        """
        Initialize a row in the technology matrix.  Each row corresponds to a reference exchange in the database, and
        thus represents a particular process generating / consuming a particular flow.

        :param flow: the LcFlow entity that belongs to the parent node's exchange
        :param process: the termination of the parent node's exchange; identity of the node. None is equivalent to a
        cutoff flow or elementary flow.  If non-null, the process must possess a reference exchange with the same flow
        or an error will be raised
        """
        self._flow = flow
        self._process = process
        if process is None:
            self._hash = (flow.get_uuid(), None)
        else:
            if len([x for x in process.reference_entity
                    if x.flow == flow]) == 0:
                raise NoMatchingReference('Flow: %s, Termination: %s' % (flow.get_uuid(), process.get_uuid()))
            self._hash = (flow.get_uuid(), process.get_uuid())

    def __eq__(self, other):
        if not isinstance(other, ProductFlow):
            return False
        return self.flow == other.flow and self.process == other.process

    def __hash__(self):
        return hash(self._hash)

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
    def is_cutoff(self):
        return self._process is None

    @property
    def is_elementary(self):
        return is_elementary(self._flow)


class ForegroundManager(object):
    """
    Class for managing a collection of linked processes.
    """
    def __init__(self, data_dir=DEFAULT_DATA_DIR):
        """

        :param data_dir:
        """
        self._data_dir = data_dir
        self._lowlink = dict()  # dict mapping product_flow to lowlink -- which is a key into TarjanStack.sccs

        self.tstack = TarjanStack()

        self._interior = []  # entries in the sparse A matrix
        self._cutoff = []  # entries in the sparse B matrix
        self._product_flows = dict()  # maps product_flow to index
        self._index = []  # maps index to product_flow

    def index(self, product_flow):
        return self._product_flows[product_flow]

    def product_flow(self, index):
        return self._index[index]

    @property
    def sccs(self):
        return self.tstack.sccs

    def scc(self, key):
        return self.tstack.scc(key)

    def lowlink(self, product_flow):
        return self._lowlink[product_flow]

    def _add_product_flow(self, pf):
        index = len(self._index)
        self._product_flows[pf] = index
        self._set_lowlink(pf, index)
        self._index.append(pf)
        self.tstack.add_to_stack(pf)

    def _set_lowlink(self, pf, lowlink):
        """
        Sets lowlink to be the lower of the existing lowlink or the supplied lowlink
        :param pf:
        :param lowlink:
        :return:
        """
        if pf in self._lowlink:
            self._lowlink[pf] = min(self.lowlink(pf), lowlink)
        else:
            self._lowlink[pf] = lowlink

    def _check_product_flow(self, flow, termination):
        """
        returns the product flow if it exists, or None if it doesn't
        :param flow:
        :param termination:
        :return:
        """
        k = ProductFlow(flow, termination)
        if k in self._product_flows:
            return self.product_flow(self.index(k))
        else:
            return None

    def _create_product_flow(self, flow, termination):
        pf = ProductFlow(flow, termination)
        self._add_product_flow(pf)
        return pf

    def add_ref_product(self, archive, flow, termination, multi_term='first', default_allocation=None):
        """
        Here we are adding a reference product - column of the A + B matrix.  The termination must be supplied.
        :param archive:
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
        j = self._check_product_flow(flow, termination)
        if j is None:
            j = self._create_product_flow(flow, termination)
            self._traverse_term_exchanges(archive, j, multi_term, default_allocation)
        return j

    def _traverse_term_exchanges(self, archive, parent, multi_term, default_allocation=None):
        """
        Implements the Tarjan traversal
        :param archive:
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
            if exch.termination is None:
                terms = [t for t in archive.terminate(exch)]
                if len(terms) == 0:
                    self.add_cutoff(parent, exch)
                    continue
                else:
                    term = ambiguous_termination(exch, terms, multi_term)
            else:
                term = archive[exch.termination]
            i = self._check_product_flow(exch.flow, term)
            if i is None:
                # not visited -- need to visit
                i = self._create_product_flow(exch.flow, term)
                self._traverse_term_exchanges(archive, i, multi_term, default_allocation)
                # carry back lowlink, if lower
                self._set_lowlink(parent, self.lowlink(i))
            elif self.tstack.check_stack(i):
                # visited and currently on stack - carry back index if lower
                self._set_lowlink(parent, self.index(i))
            else:
                # visited, not on stack- nothing to do
                pass
            self.add_interior(parent, exch, i)

        # name an SCC if we've found one
        if self.lowlink(parent) == self.index(parent):
            pass
        self.tstack.label_scc(parent.key)

    def add_cutoff(self, parent, exchange):
        """
        Create an exchange for a cutoff flow (incl. elementary flows)
        :param parent:
        :param exchange: an LcExchange belonging to the parent node
        """
        self._cutoff.append(MatrixEntry(parent, None, exchange))

    def add_interior(self, parent, exchange, term):
        """

        :param parent:
        :param exchange:
        :param term:
        :return:
        """
        self._interior.append(MatrixEntry(parent, term, exchange))
