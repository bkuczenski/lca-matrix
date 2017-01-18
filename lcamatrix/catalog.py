"""
A collection of helper functions for extracting + loading catalogs.  Some things that are required include:

 * specifying characterization factors
 * constructing market / mixer processes
 * applying allocation

Other things that may also be required, but that don't yet have a clear mechanism:

 * identifying synonyms
"""


from collections import namedtuple
from lcatools.entities import MissingFactor

'''
A ConfigFlowCharacterization provides a procedural mechanism for specifying flow quantity characterizations after
loading an archive.  The 'flow_ref' and 'quantity_ref' have to lookup successfully in the archive.  Note- provisioning
the archive itself is ad hoc. Don't have a structural way to do that yet, but someday it could be done with linked
data (maybe?)
'''
ConfigFlowCharacterization = namedtuple("ConfigFlowCharacterization", ('flow_ref', 'quantity_ref', 'value'))

'''
A ConfigAllocation provides a procedural mechanism for specifying quantity-wise allocations of processes at load time.
All that is required is a quantity; the process knows how to perform the allocation.  Note that reference flows not
characterized with respect to the quantity will receive zero allocation.  So apply_flow_config first.
'''
ConfigAllocation = namedtuple("ConfigAllocation", ('process_ref', 'quantity_ref'))


def apply_flow_config(archive, flow_characterizations, overwrite=False):
    """
    Applies a list of ConfigFlowCharacterizations to an archive.
    :param archive:
    :param flow_characterizations:
    :param overwrite: [False] overwrite if characterization is present
    :return: nothing- the changes are made in-place
    """
    for cfc in flow_characterizations:
        if not isinstance(cfc, ConfigFlowCharacterization):
            raise TypeError('Entry is not a ConfigFlowCharacterization\n%s' % cfc)
        flow = archive[cfc.flow_ref]
        qty = archive[cfc.quantity_ref]
        if flow.has_characterization(qty):
            if overwrite:
                flow.del_characterization(qty)
            else:
                print('Flow %s already characterized for %s. Skipping.' % (flow, qty))
                pass
        flow.add_characterization(qty, value=cfc.value)


def apply_allocation(archive, allocations, overwrite=False):
    """
    Applies a list of ConfigAllocations to an archive.

    If overwrite is True, the process's allocations are first removed.

    If overwrite is False, the process is tested for allocation under each of its reference flows- if any are
    already allocated, the allocation is aborted for the process.

    :param archive:
    :param allocations:
    :param overwrite: [False] whether to strike and re-apply allocations if they already exist.
    :return:
    """
    for al in allocations:
        if not isinstance(al, ConfigAllocation):
            raise TypeError('Entry is not a ConfigAllocation\n%s' % al)
        p = archive[al.process_ref]
        qty = archive[al.quantity_ref]
        is_alloc = False
        if overwrite:
            for rf in p.reference_entity:
                p.remove_allocation(rf)
        else:
            for rf in p.reference_entity:
                try:
                    is_alloc |= p.is_allocated(rf)
                except MissingFactor:
                    is_alloc = True
                    break

        # now apply the allocation
        if is_alloc:
            print('Allocation already detected for %s. Skipping this configuration.' % p)
            continue
        else:
            p.allocate_by_quantity(qty)
