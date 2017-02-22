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

        self._hash = (flow.uuid, None)
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
        """
        Product flow key is (uuid of reference flow, uuid of process)
        :return:
        """
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

    def __str__(self):
        return '%s==%s' % (self._process, self._flow)
