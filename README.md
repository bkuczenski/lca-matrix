# LCA Matrix Work

Tasks:

 * replicate Partial Ordering work:
   1. Construct a square technology matrix from a collection of inventories BA^-1 * x
   2. Modify the matrix to be a normalized adjacency matrix B(I-A)^-1 * x
   3. Detect strongly connected components
   4. Order the matrix
   
 * Implement foreground / background discovery
 
 * Construct foreground studies from final demand
 
 * Construct case studies for publishing paper.
 
 
Ideally this will be done without using `lca-tools`, but on the other hand, why should that be an ideal?

Upon review, given that the linking algorithm must be capable of auto-creating new market processes, it only makes sense for `lca-tools` to be an import.  That is ultimately a good thing, because it will provide guidance in instrumenting `lca-tools` with unit tests.

Other tasks that are blocking me:

 - lack of testing / interface definition for `lca-tools`, particularly with exchanges, especially with allocated exchanges
 - Compartments refactor needs to be finished and merged.
 - flow synonym problems with USLCI-- looks like they don't exist! (at least when we limit our concern to interior flows-- see narrative below) 
 - `lca-tools` should not depend on `pandas` - I'm only using it in one place, and only for spreadsheet handling.

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
  
  * Processes with multiple reference flows, all of which appear in the interior matrix, will require allocation during linking.
  
  * Processes with multiple reference flows, of which at least one does not appear in the interior matrix, will be allocated over the interior flows only.  **Problem** (for which the solution is to add a market process
  


#### Selecting among multiple producers

In ecospold v2 (currently exclusively used to its full technical extent by ecoinvent), each interior flow features at least one _market_ process, which is a mixer of various production sources to form a marketed commodity for a specific region.  Market processes are an organizational tool that assigns production mix and transport requirements to every commodity.  Ecoinvent also uses "market groups" which are further idiosyncrasies of the ecospold-ecoinvent system, and which organize subsidiary markets into larger-scoped aggregations.

The existence of market processes is ensured through the use of linking software such as [ocelot](http://ocelot.space) (currently also exclusively used for ecoinvent).  However, absent the expansion of scope of ocelot to include non-ecoinvent systems, this will have to be handled manually.




I should be able to build a minimal technology matrix using _only_ interior flows!

then I can move forward on Tarjan without resolving the compartments refactor.

Here's how I proceed:

For USLCI, with strong process-reference flow identity:

 * Identify interior flows -- create a list, in arbitrary order
 * generate a 
 
 
