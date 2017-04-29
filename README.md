# LCA Matrix Work

This repository provides software for matrix-based computation of LCI results.

## How does it work?

The BackgroundManager software performs a [partial ordering](https://www.researchgate.net/publication/282879534) of a database stored in an [LcArchive](https://github.com/bkuczenski/lca-tools).  The ordering identifies strongly connected components, and thus allows for automatic distinction between the foreground and background of an LCI database:

 * the _background_ includes processes whose parameters affect all results in the database;
 * the _foreground_ includes processes whose parameters affect only the product system that uses them.

The ordering results in an implicit grouping of the _flows_ in a database:

 * _product flows_ are the reference flows of a process in the database.
 * _exterior_ flows are flows for which only one terminus is present in the database (e.g. flows that appear as an input but nowhere as an output, or vice versa).  Reference flows are excluded from exterior flows.
   * Exterior flows whose _compartments_ are elementary are called elementary flows or _emissions_ (even if they are inputs).
   * Intermediate exterior flows are called _cutoff flows_.

Mathematically, product flows make up the rows and columns of the $A$ matrix, and exterior flows make up the rows of the $B$ matrix.  Although the $B$ matrix is traditionally used only for emissions, here they are included in the $B$ matrix together because emissions and cutoff flows are mathematically identical.

### Model components

The `BackgroundManager` satisfies the LciInterface specified in the `lca-tools` repository (TODO).  


### Allocation Issues

Conventional practice in LCA requires multi-functional processes to be represented as single-output processes in order to construct an invertible technology matrix.  However, there are alternative approaches that still allow an invertible matrix to be constructed.

 1. In-place substitution.  If a co-product is substituted by a single-output process that already exists in the database, then the two products can simply be made synonymous -- both assigned the same row+column in the technology matrix.  In this case, inversion of the matrix can compute the activity levels of all processes in a way that solves the production problem, with the single-output process making up the difference (positive or negative) between the joint production amount and the demanded amount.  This is (or appears to be effectively) what Ecoinvent's Allocation-at-the-point-of-substitution system model implements.

 2- Surplus Co-product.  An un-allocated multi-output process can be implemented directly as-is in a technology matrix simply by selecting any one of the reference exchanges to be the process's single "primary" reference, and listing the process's other reference flows directly (negating outputs as appropriate).  The technology matrix will still be invertible if the other reference flows are also provided "cut-off" processes to act as unconstrained sources or sinks for the co-products.  In this case, the unallocated process will simply be operated at the level required by demand for its primary reference product, and the other reference products will be produced as surplus products.  Demand for these products from elsewhere in the system will reduce the surplus, and the matrix inversion will end up computing the net surplus (or net requirement, if negative) of each surplus co-product.  Co-products in this arrangement do not balance when the technology matrix is inverted-- but that may be appropriate if the consuming processes are not modeled within the database.  

Implementing co-products-to-cutoff requires knowledge of which materials are surplus products. This cannot be determined in advance of constructing the entire technology matrix, and it may vary depending on the functional unit required.  Currently, the default behavior is to apply the surplus-coproduct approach when unallocated processes are encountered, but [TODO: keeps a record of surplus co-products, which are accounted for in the exterior matrix.]  These will always be a subset of cutoff flows.  [TODO: The BackgroundManager also has a flag `fully_allocated` which will report False if the list of surplus co-products is nonempty.]


## narrative

First, thinking about interior flows / terminations / linking:

USLCI has a strong identity matching process to reference product (each product is uniquely produced by one process)

Actually, let's test that

### USLCI count MAKE table redundancy

On review, there are a number of elementary flows that appear as both inputs and outputs- these must necessarily be excluded from the interior flows.


Let's check to see if the USLCI data from ecospold has the same synonym problems I encountered before.

Looks like NO!  In particular, there are no duplicates containing the terms 'diesel', 'gasoline', 'refinery', or 'electricity' (on case-insensitive search, via `lower()`).

That could make things a bit easier.

## USLCI load-in




### USLCI Synonym script

This is a cheesy point-test, and not a rigorous effort to detect all possible synonyms.

	from lcamatrix import from_json, make_flows, interior_flows
	
	CATALOG = '/data/GitHub/lca-catalog/catalogs/'
	USLCI = 'uslci_ecospold.json.gz'
	
	# build catalog
	us = from_json(os.path.join(CATALOG, USLCI))

	# us_flows is a dict of entityId -> flow dict
	us_flows = make_flows(us)

	# in_us is a set of flows (by id) that show up as both input and output exchanges
	in_us = interior_flows(us)

	# elementary search
	sorted([us_flows[f]['Name'] for f in in_us if us_flows[f]['Name'].lower().find('diesel') >= 0])
	
	# inspect the outcome to find repetitions
	
	
### Building a technology matrix

#### Identifying commodities

I can divide the flows up by type:

 1. Use a standalone Compartment manager to distinguish between product and elementary flows
 2. Within product flows, detect _interior_ flows-- which appear as both an input and an output.
 3. Non-interior product flows are used as either an output or an input, but not both
   a. exterior inputs are _cutoff_ flows
   b. exterior outputs are _reference products_

(So already we are blocking on the compartments refactor.)

The interior flows are the only things required to be included in the technology matrix-- cutoffs and reference products can be added on later.  We will refer to this scope-reduced technology matrix as the _interior matrix_.

Once the interior flows are identified, it becomes possible to construct a list of _commodities_ which will be the rows and columns of the interior matrix.  The rows and columns of the interior matrix are determined as follows:

 * for each process:
   1. Identify the process's _reference exchanges_.
   2. If the process has no designated reference exchanges, take the intersection of its outputs with the set of interior flows.  (This assumption bears revisiting later on).
   3. Group the commodities by interior flow (in the currently imagined python implementation, interior flows are organizational objects and so they group themselves.  The process would register itself with the interior flow as a producer).  Note that the direction of flow is _not_ assumed here-- reference flows can be either inputs or outputs, as long as the same interior flow is not a reference input and a reference output in the same database.
   4. Each unique combination of a process and an interior flow becomes a commodity in the matrix.
   5. (conjectural) each reference exchange found not to be interior should be added to the interior flows and registered with the process.

Interior flows can now be used to construct an interior matrix by matching each process's _non-reference_ exchanges to a complementary exchange known to an interior flow.  In the case of a linked ecospold v2 system, we are done-- the exchanges can now be transformed into an adjacency matrix.

However, in the case of ILCD or ecospold v1 formatted data, there are [at least?] two possible confounding factors that challenge the creation of automatically linked system models:

 * __multiple reference flows__: For a process that generates multiple reference flows, it is necessary to determine how to allocate the process's dependencies and emissions among the products.  Ecoinvent system models implement various allocation strategies, so this is only a challenge when dealing with other data sources.

 * __multiple producers__: In the event that several production models all produce the same product, it is insufficient to only specify the interior flow being exchanged- the dependent process must also specify the producing process.  Ecospold v2 exchanges include `activityLinkId` fields which specify the partner process; but no other formats support this.


#### Multiple reference flows

The co-product allocation problem is perhaps the most persistent methodological challenge in life cycle assessment.  It has been conclusively determined that no automatic, procedural allocation scheme is appropriate _in general_ to all cases (e.g. Pelletier 2015 for a comprehensive review).

In this application, the allocation challenge is reduced _somewhat_ because the interior matrix only includes interior flows, which are flows that appear as both inputs and outputs in the database.

  * Processes with one or more reference flows, none of which appear in the interior matrix, will not themselves be included in the interior matrix and will only exist as foreground processes.  In this case, the data user can handle allocation at the time of use.
  
  * Processes with one reference flow, which appears in the interior matrix, will not require allocation.
  
  * Processes with multiple reference flows, all expected to appear in the interior matrix, will require allocation during linking.
    
I will need to be able to specify allocations later on.

#### Multiple producers

In ecospold v2 (currently exclusively used to its full technical extent by ecoinvent), each interior flow features at least one _market_ process, which is a mixer of various production sources to form a marketed commodity for a specific region.  Market processes are an organizational tool that assigns production mix and transport requirements to every commodity.  Ecoinvent also uses "market groups" which are further idiosyncrasies of the ecospold-ecoinvent system, and which organize subsidiary markets into larger-scoped aggregations.

The existence of market processes is ensured through the use of linking software such as [ocelot](http://ocelot.space) (currently also exclusively used for ecoinvent).

In our case we can either assume a def




I should be able to build a minimal technology matrix using _only_ interior flows!

then I can move forward on Tarjan without resolving the compartments refactor.

Here's how I proceed:

For USLCI, with strong process-reference flow identity:

 * Identify interior flows -- create a list, in arbitrary order
 * generate a 
 
 
