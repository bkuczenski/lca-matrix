from collections import defaultdict

from lcamatrix.product_flow import ProductFlow


class TarjanStack(object):
    """
    Stores the current stack and provides a record of named SCCs
    """
    def __init__(self):
        self._stack = []
        self._stack_hash = set()
        self._sccs = defaultdict(set)  # dict mapping lowest index (lowlink = SCC ID) to the set of scc peers
        self._scc_of = dict()  # dict mapping product flow to SCC ID (reverse mapping of _sccs)

        self._component_cols_by_row = defaultdict(set)  # nonzero columns in given row (upstream dependents)
        self._component_rows_by_col = defaultdict(set)  # nonzero rows in given column (downstream dependencies)

        self._background = None  # single scc_id representing largest scc
        self._downstream = set()  # sccs on which background depends
        self._fg_index = dict()  # maps product_flow.index to a* / b* column -- STATIC
        self._bg_index = dict()  # maps product_flow.index to af / ad/ bf column -- VOLATILE

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
        if ml > 1:
            self._background = ind
            self._set_downstream()
            self._generate_bg_index()

    def _set_downstream(self, upstream=None):
        """
        recursive function to tag all nodes downstream of the named node.
        :param upstream: [None] if none, use background
        :return:
        """
        if upstream is None:
            if self._background is None:
                return
            upstream = self._background

        for dep in self._component_rows_by_col[upstream]:
            if dep != upstream:  # skip self-dependencies
                self._downstream.add(dep)
                self._set_downstream(dep)

    def _generate_bg_index(self):
        if self._background is None:
            self._bg_index = dict()
        else:
            bg = [b.index for b in self.background_flows()]  # list of ProductFlow.index values of background nodes
            self._bg_index = dict((ind, n) for n, ind in enumerate(bg))  # mapping of *pf* index to a-matrix index

    def _generate_foreground_index(self):
        """
        Perform topological sort of fg nodes
        :return:
        """
        fg_nodes = set()
        fg_ordering = []
        for k in self._sccs.keys():
            if k != self._background and k not in self._downstream:
                if len(self._component_cols_by_row[k]) == 0:  # no columns depend on row: fg outputs
                    fg_ordering.append(k)
                else:
                    fg_nodes.add(k)

        while len(fg_nodes) > 0:
            new_outputs = set()
            for k in fg_nodes:
                if len([j for j in self._component_cols_by_row[k] if j not in fg_ordering]) == 0:  # all upstream is fg
                    new_outputs.add(k)
            fg_ordering.extend(list(new_outputs))  # add new outputs to ordering
            for k in new_outputs:
                fg_nodes.remove(k)  # remove new outputs from consideration

        self._fg_index = dict((ind, n) for n, ind in enumerate(fg_ordering))

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
        # self._generate_foreground_index()

    @property
    def background(self):
        return self._background

    @property
    def ndim(self):
        return len(self._bg_index)

    @property
    def pdim(self):
        return len(self._fg_index)

    def foreground(self, index):
        """
        returns a list of foreground SCCs that are downstream of the supplied scc ID (itself included)
        This is a list of components necessary to build the foreground fragment tree
        :param index:
        :return: topologically-ordered, loop-detecting list of non-background SCC IDs
        """
        if self.is_background(index):
            return []
        queue = [index]
        fg = []
        while len(queue) > 0:
            current = queue.pop(0)
            if current not in fg:
                queue.extend([k for k in self._component_rows_by_col[current] if not self.is_background(k)])
                fg.append(current)
        return fg

    def foreground_flows(self, outputs=False):
        """
        Generator. Yields product flows in
        :param outputs: [False] (bool) if True, only report strict outputs (nodes on which no other nodes depend)
        :return:
        """
        for k in self._sccs.keys():
            if k != self._background and k not in self._downstream:
                if outputs:
                    if len(self._component_cols_by_row[k]) > 0:
                        continue
                for pf in self.scc(k):
                    yield pf

    def background_flows(self):
        """
        Generator. Yields product flows in the db background or downstream.
        :return:
        """
        if self._background is None:
            return
        for i in self.scc(self._background):
            yield i
        for i in self._downstream:
            for j in self.scc(i):
                yield j

    def bg_dict(self, pf_index):
        """
        Maps a ProductFlow.index to the row/column number in A* or column in B*
        :param pf_index:
        :return:
        """
        try:
            return self._bg_index[pf_index]
        except KeyError:
            return None

    def fg_dict(self, pf_index):
        """
        Maps a ProductFlow.index to the column number in Af or Ad or Bf
        :param pf_index:
        :return:
        """
        try:
            return self._fg_index[pf_index]
        except KeyError:
            return None

    def is_background(self, index):
        """
        Tells whether a Product Flow index is background. Note: SCC IDs are indexes of the first product flow
        encountered in a given SCC, or the only product flow for a singleton (i.e. acyclic foreground) SCC
        :param index: product_flow.index
        :return: bool
        """
        return index in self._bg_index

    def scc_id(self, pf):
        return self._scc_of[pf]

    def sccs(self):
        return self._sccs.keys()

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
