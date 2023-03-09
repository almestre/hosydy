# hosydy

## General description

Model of adaptive colonisation across a parasitism–mutualism gradient.

This code was created to implement an individual-based model to investigate adaptive colonisation by hosts and their symbionts across a parasite–mutualist continuum. The host must adapt in order to establish in the novel habitat, and the symbiont must adapt to track evolutionary change in the host. We introduce variation in dispersal at two scales: host migration and symbiont transmissibility.

The code was used to perform the analyses included in a MS submitted to a journal:
Mestre A, Butlin RK, Hortal J and Rafajlović M (2023). Adaptive colonisation across a parasitism-mutualism gradient. Submitted.

## Details of the code that implements the model

For this paper, we wrote an individual-based model using C++ language (available in the repository …). We highlight here four main features of the code. First, it is based on object-oriented programming techniques suitable to reflect the hierarchical nature of our system. Second, individuals and populations are stored in fixed-size vectors that optimise execution times. Third, ‘slot map’ data structures provide fast and secure access to individuals by unique identifiers, thus optimising communication among objects —e.g., a symbiont accessing to information from its host to calculate its fitness, or communication between a newborn host and its parents for implementation of vertical transmission. Fourth, genetic information is stored in bits and managed through bitwise operations, thereby optimising memory usage. Next, we provide some details of each feature.

### Object-oriented programming

The code that implements our model follows a nested structure with the following classes:

#### Gamete class
This class instantiates objects representing individual gametes with their haploid genotype.

#### Host gamete class
A specialisation of the Gamete class for host gametes. It also stores the ID of the original host individual that produced the gamete, so that algorithms of vertical transmission can locate parents using the information stored in both gametes that will lead to a newborn.

#### Organism class 
This class stores and manages sex, genotype and phenotype of an individual.
It contains the algorithms that create a new haploid genotype from the diploid genotype by free recombination (for gamete production).

#### Host class
Inherited from Organism class (inherits both data and functionality)
This class instantiates objects representing individual hosts.
Each individual host stores:
    • Sex
    • Genotype
    • Phenotype
    • Host ID, a 64-bit number acting as a unique identifier
    • Symbiont local population size, i.e., the number of symbionts it harbours
    • An object representing a population of symbionts (see below).
Moreover, all the objects of the Host class share a unique copy of the following data (which are assumed to be the same for all the individuals):
    • Carrying capacity of a local population of symbionts (kL)
    • Maximum intrinsic population growth rate (Rmax)
    • Strength of selection (Vs)
    • Growth rate of the interaction effect size (b)
    • Host death rate per capita per reproduction event (μH).
In addition, it contains the associated functionality to manage these data, including gamete production depending on its current state (i.e., fitness).

#### Symbiont class
Inherited from Organism class (inherits both data and functionality)
This class instantiates objects representing individual symbionts.
Each individual symbiont stores:
    • Sex
    • Genotype
    • Phenotype.
All the objects of the Symbiont class share a unique copy of the following data (which are assumed to be the same for all the individuals):
    • Maximum intrinsic population growth rate (rmax)
    • Strength of selection (Vs)
    • Per capita per generation emigration rate for vertical transmission (γ).
In addition, it contains the associated functionality to manage these data, including gamete production depending on its state (i.e., fitness).

#### Population class template
This class instantiates objects that represent a population of either hosts or symbionts.
Each population-type object contains the following data: 
    • A vector of individuals (either hosts or symbionts)
    • Patch ID, either ID of the land patch inhabited by a host population (for hosts), or the ID of the host inhabited by a symbiont infrapopulation (for symbionts)
    • The current (infra)population size
    • Population ID.
Population-type objects have available all the functionality for implementing population-level processes including immigration, reproduction and mortality, and calculating the key output variables involved in these processes.
In the case of host populations, because we need access to hosts by ID, the class template variant for hosts includes a ‘slot map’ that manages the vector of host individuals, and the associated functionality (see next section).
By contrast, the class template specialisation for symbionts includes the ID of the host harbouring the population symbionts (for information transfer purposes).

#### Metapopulation class template
This class instantiates objects that manage a vector of either symbiont or host populations. For our current research question, we only use the symbiont-type template specialisation, which instantiates objects representing a global population of symbionts (creating host metapopulations is also possible with our code but this option is not utilised here). An object of this type stores a vector of symbiont infrapopulations, and the associated functionality for implementing processes acting at the symbiont global population level, including reproduction, vertical transmission, creation of new infrapopulations by host immigration from the continent, or destruction of infrapopulations by host mortality events; it also includes the functionality for calculating the key output variables involved in these processes.
As mentioned above, each symbiont population stores the ID of its host, which is unique for each symbiont population at a given time step. Thus, a symbiont-metapopulation object instantiated by this class stores a slot map that manages fast access to symbiont populations based on their host’s ID, thus enabling the algorithms of this class to implement complex processes that require information transfer among objects of different types (e.g., vertical transmission).

#### Patch class
This class instantiates objects representing a land patch that may serve as a living place for a population of hosts.
Each land patch contains:
    • Patch carrying capacity for hosts
    • The optimal phenotype for hosts living in the patch
    • A patch ID
    • The ID of the host population inhabiting the patch.
Objects of patch class have the functionality required to manage their data, provided by the class. In this study we only implement a situation with a single patch representing the island.

#### Rng class
This class contains all the functionality used by the algorithms of other classes to generate random numbers based on a variety of distributions, including uniform, Bernoulli, binomial, Poisson, normal and negative binomial.
The class instantiates and stores a pseudo-random number generator that produces 32-bit pseudo-random numbers using the Mersenne twister algorithm (https://cplusplus.com/reference/random/mt19937/).
It also manages input and output of DAT files that store the states of the pseudo-random engine, so that replication of simulations obtaining identical results is possible.

#### Parameter class
This class manages parameter setting based on both default values and an input file of JSON type.

#### Output class
This class manages the creation of the output files that will store the simulation data for subsequent analyses. 

#### Simulation class
This class manages the parameters and functionality that controls the simulation procedure.
The class stores:
    • Number of symbiont cycles (i.e., steps) per host reproductive event
    • Number of host reproductive events per year
    • Number of years simulated
    • ID of the simulation scenario
    • ID of the simulation replicate
    • The algorithm that runs the simulations.

### Fixed-size vectors to store individuals and populations

In our code, population- and metapopulation-type objects use fixed-size vectors to store individuals and populations. We take advantage of the fact that we have an idea about the population-size limits provided by carrying capacities (i.e., K + some additional amount). In this way, a vector of individuals stored within a population-type object is initialised with a maximum potential number of individuals estimated from carrying capacities (which can be tested by preliminary simulations). Then, we simulate variation in population size without the need for object constructions/destructions (which is time-consuming). We do that by modifying the information stored within the objects without destroying them (including swaps and “reinitialisations”). The basic idea is that we “recycle” objects, and only the first N objects of the vector are being "used" at a given time step (where N is the population size). When a new individual is born, the corresponding algorithm reinitialises the object at position N (vector indices in C++ start from 0); and increments the variable N (stored in the population-type object) by one. When an individual at given position dies, the corresponding algorithm swaps the object in that position with the object at position N-1 (i.e., the last object that is being used); and decrements N by one (the same applies to population-type vectors stored within metapopulation-type objects).

Using fixed-vector sizes is not optimal in terms of memory usage, but it is faster (resizing vectors of objects is costly). So, here we decided to optimise execution time rather than memory usage (space-time trade-off). The approach required a method to manage IDs of individual/population vectors, explained in the next section.

### A “Slot map” structure to manage IDs

Many algorithms in our model required a method for searching an object by its ID across a vector that stores objects of its type (e.g., searching for a given individual host stored in a host-population vector). However, we wanted to avoid loops through long vectors to find objects by ID, or using higher-performance algorithms that require prior sorting of elements. We found an optimal and quite simple solution from the field of game design: a "slot map" data structure (https://seanmiddleditch.github.io/2013-01-05-data-structures-for-game-developers-the-slot-map/)

A slot map has three components:
a) A vector with the objects
b) An indirection list
c) A free list.

The indirection list is a vector of integers that contains the information for fast access to objects based on their ID. The free list is a stack (last-in-first-out data container) that stores the indices that are not being used.

We want IDs to act as recyclable indices, but at the same time it is worth having unique identifiers that are not recycled. The solution is to split the ID into two parts, an index (which is recycled) and a version, so that the ID as a whole is unique and never recycled. Objects store the whole ID, whereas both the indirection and free list work based on the index part. For instance, let’s suppose we have a vector storing a population of hosts that is managed by an indirection list and a free list, as follows:

Vector of objects (Population):
obj0 { ID= (0, 0) } // ID = (a, b) where a = index and b = version
obj1 { ID= (0, 0) }
obj2 { ID= (0, 0) }
obj3 { ID= (0, 0) }
obj4 { ID= (0, 0) }
obj5 { ID= (0, 0) }

Indirection list (IndirectionList): 0 1 2 3 4 5

Free list: 5 4 3 2 1 0 
// elements here are always taken from and added to the back (i.e., stack container)

N = 0

In our example we have 6 objects. The N value indicates that the population size is 0, so that all the objects are unused. After the birth of four hosts, the data would be as follows:
 
Vector of objects (Population):
obj0 { ID= (0, 1) } // ID = (a, b) where a = index and b = version
obj1 { ID= (1, 1) }
obj2 { ID= (2, 1) }
obj3 { ID= (3, 1) }
obj4 { ID= (0, 0) } // object never used
obj5 { ID= (0, 0) } // object never used

Indirection list (IndirectionList): 0 1 2 3 4 5

Free list: 5 4

N = 4

The free list indicates which is the next "available" ID_index to use (here the 4). Notice that each time an object is used, the following actions are performed:
    • Its version is incremented by 1
    • Its ID_index is removed from the free list
    • N is incremented by 1.
After that, N indicates that we have 4 hosts, which we know are in positions 0-3 of the objects' vector.

Let’s suppose that the host with ID = (1, 1) dies. This host is stored in Obj1. The required actions are:
    1. Swap Obj1 with Obj3 (i.e., the object at position N-1)
    2. Swap the corresponding values in the indirection list
    3. Increment ID_version of Obj1
    4. Add the ID_index to the free list for its reuse
    5. Decrement N by one

The resulting data are:

Vector of objects (Population):
obj0 { ID= (0, 1) } // ID = (a, b) where a = index and b = version
obj3 { ID= (3, 1) }
obj2 { ID= (2, 1) }
obj1 { ID= (1, 2) } // currently unused
obj4 { ID= (0, 0) } // object never used
obj5 { ID= (0, 0) } // object never used

Indirection list (IndirectionList):
0 3 2 1 4 5

Free list:
5 4 1

N = 3

As it can be seen, we have now three living hosts. The obj1 (with ID_index = 1) is currently unused. But, because the free list is a stack container (last-in-first-out), it will be the next object to be "recycled" in the next birth. Also, note that the version-part of the IDs informs about whether the object is currently being used (uneven number) or unused (even number), and how many times it has been used.

Let’s suppose a given algorithm requires access to a host with a given ID. This search is implemented by a “searching” algorithm that uses the index part of the ID and the indirection list to access the host as follows:
ObjectLocation = IndirectionList[ID_index]
Objects[ObjectLocation]

In our example, the algorithm would get access to the host with ID = (3, 1) as follows:
ObjectLocation = IndirectionList[3] // which equals 1
Objects[1] // i.e., the vector address for obj3, which has the ID = (3, 1)

Then, the algorithm uses the whole ID (index + version) to make sure that the reached object stores the correct expected host before making actions:
if (ID_stored_in_the_located_object == ID_stored_in_the_symbiont) {do whatever you wanted}

After locating the target object, the “searching” algorithm creates a pointer to the object, and passes the pointer to other algorithms that use the pointer to access the object and destroy the pointer following usage. This is a fast and secure method to access individuals by ID.

### Genetic information stored in bits

The memory costs of storing genetic information for hosts and symbionts can be a concern. Our initial idea was to use a string to store the genotype. A string is a vector of characters, where each character can represent an allele. The size of a character type is equal to 1 byte (i.e., 8 bits). Thus, the size of a string-format genotype is at least 2*L bytes (i.e., 2*L*8 bits; where L is the number of biallelic loci). The size of a boolean type is also 1 byte (thus being equivalent to character type in terms of memory usage). Therefore, to save memory, we decided to store allellic information in bits. 

Our code interprets bits as boolean variables, and uses bitwise operations (bit shifts and masks) for direct manipulation of the bits that code for a variable. The code uses a 64-bit integer to store a genotype with L = 32. A 64-bit integer can be split into two parts:
    1. A lower 32-bit integer
    2. An upper 32-bit integer
where each part represents one of the two alleles of the L loci. Thus, we use 64 bits (i.e., 8 bytes) to store a genotype. The string-format would require 2*32*8= 512 bits (i.e., 64 bytes). In this way, we get a 1/8 reduction in memory usage. Moreover, manipulations of objects (e.g., initialisation or swapping) that contain an integer-format genotype is faster, compared with those containing a string-format genotype.

In fact, our code also represents unique IDs by 64-bit integers (see previous section). In this case, the lower 32-bit part is the index and the upper-32-bit part is the version. Then, bitwise operations are implemented to replace the index part and increment the version. The index alone can also be represented by a 32-bit integer —which is the format of the elements stored in the indirection and free lists.
