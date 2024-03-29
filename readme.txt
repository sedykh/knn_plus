knn+ is a k-nearest neighbor modeling method with descriptor selection. 
Its input is a "split" file with a list of training and test sets (as generated by datasplit)

Example of classification modeling: 
'knn+ demo_input_split -M=CTG -O=SA -LOGALL' 
(check generated log file for the progress)


Modeling output is an ensemble of individual kNN models (.tbl file), 
each defined by its subspace (set of descriptors) and k parmeter.
The file can be parsed or used entirely for the consensus prediction.

----------------------------------------------------------------
kNN+ V3.0 - Integrated family of kNN algorithms. 2014
Usage: 'knn+ inputfile [flags]'
Allowed input files: .t2t, .list (for modeling)
 and .tbl, .list (for prediction)

Use '-4PRED=testfile' for prediction mode, otherwise it is model building mode.
'-OUT=...' - user-defined output file, otherwise inputname is used

'-O=...' - descriptor selection methods: 'SA'<def.>,'GA','GA1','AC','PS','++','--','XX'
'GA','GA1' - genetic algorithms, 'SA' - simulated annealing
'AC' - ant colony, 'PS' - particle swarm
'++','--' - stepwise kNN, 'XX' - combinatorial search
For more info type 'knn+ -O=...' with the method of interest.

'-D=...' - range of dimensions to use, e.g. '-D=5@50' or
'-D=5@50@2' to scan with steps of 2 <def.>

'-KR=...' sets knn range (#nearest neighbors) or R(Radial cut-off)
e.g. '-KR=1@9' <def.>(uses 1-9 neighbors) or '-KR=0.5'

'-M=...' - model type: 'CNT' continuous <def.>,'CTG' category,'CLS' classes

'-WT=...' - choosing and weighting neighbors; e.g. -WT=AE2.0
'2.0' is power-coefficient of weight-f(); E - exponential weight-f()
'H' - hyperbolic, 'M' - Minkowski-kind weight-f() <def.>,
'N' - no weight-f(), all neighbors are equal;
'B' - use absolute distances in weight-f(), 'R' - relative <def.>,
'S' - use direct-distances, 'Q' - squared distances <def.>;
'A' - include all qualified neighbors, 'K'- strictly k only,
'V' - extra neighbors share votes <def.>,
'F' - check extra neighbors in full-D space

'-LGO=' - size or fraction of Leave-Group-Out self-prediction <def. 1>
'-EVL=' - models filter by train/test results: -EVL=0.5@0.6 <def.>
For continuous kNN it means q2 >0.5 and R2>0.6
A - alternative fit-index; E - error-based fit-index
V - aver.error based (only for discrete activity modeling)
R - use R2 as a check for test set <def.>
F - use external-Q2 for test set
Q - use the same index for training and test sets
S - skips post-evaluation (no check of correlation slope, etc.).
For more info type 'knn+ -EVL='.

'-AD=' - applicability domain: e.g. -AD=0.5, -AD=0.5d1_mxk
'0.5' is z-cutoff <def.>; d1 - direct-distance based AD <def. is dist^2>
Additional options of AD-checking before making prediction:
'_avd' - av.dist to k neighbors should be within AD (traditional)
'_mxk' - all k neighbors should be within AD
'_avk' - k/2 neighbors within AD, '_mnk' - at least 1 within AD <def.>

'-Z=zx' - metric f() to use; x - is a power coff. <def. 2.0>
z: E -Euclidean <def.>, T -Tanimoto, R -Corr., C -Cosine

'-SEED=' to start from given dimensions (e.g. =dscr1@dscr2@... or =filename.mod)

'-SRND=..': to seed randomizer <default is by time>
Specific settings for categories/classes kNN:
'-N=' - sets number of groups and their weights and break-points
e.g. '-N=3W0.6@0.2@0.2B0.5@1.5' sets 3 groups with 0.5 and 1.5 breakpoints and some weights
'-PNL=' - penalty factor for unbalanced prediction <def. 0>
'-SEP=' - weight of class-separation term in model evaluation <def. 0>

'-LAXDIMS' - predicted file's descriptors can be in any order, but labels must match
'-2OLD' - save individual models into .mod files (old format)
'-LOGALL' - report progress of each optimization cycle
'-DOMODS' - report comparison of individual models
'-DONEIBS' - report neighbors for each prediction

-------------------------------------

---Simulated Annealing settings: (knn+ -O=SA)
'-SA@...' - Simulated Annealing settings: e.g. '-SA@B=3@TE=-2@K=0.6@DT=-3@ET=-5'
'..@N=' - #SA runs to repeat; '..@D=' - #mutations at each T
'..@FULL' - to do all mutations at each T; '..@B=' - #best models to store
'..@T0=x' - start T (10^x); '..@TE=' - final T; '..@DT=' - convergence range of T
'..@K=' - T decreasing coeff.; '..@M=' - mutation probability per dimension



---Genetic Algorithm settings: (knn+ -O=GA or GA1)
GA: #genes in solution is #descriptors, each gene is 0 or 1.
GA1: gene values are dim ids, solution size is the allowed range of dimensions.

'-GA@...' - Genetic Algorithm settings: e.g. '-GA@N=500@D=1000@S=20@V=-4@G=7'
'..@N=' - population size; '..@D=' - max.#generations; '..@I=' - ideal fit
'..@S=' - #stable generations to stop; '..@V=' - min signif fitness diff (10^x)
'..@X=' - crossover rate; '..@M=' - mutation rate
'..@C=' - crossover mode: 1P, 2P, UN, 12
(i.e. ONE_POINT, TWO_POINT, UNIFORM, TWO_OR_ONE_POINT modes)
'..@P=' - parent selection mode: RANK, TOUR, RLTT
'..@G=' - group size for tournament ('TOUR') selection of parents
'..@E=' - to retain best solutions; e.g. '@E=0.01' (population portion)
or '@E=7' (#solutions). Use '@E=OFF' to disable, default is ON with 0.01
'..@Z=' - to set a penalty term for the solution size <default is 0.1>



---Ant Colony settings: (knn+ -O=AC)
'-AC@...' - Ant Colony settings: e.g. '-AC@N=300@D=600@S=20@V=-4@E=3@P=0.1'
'..@N=' - population size <def.100>; '..@D=' - max.#cycles <def.300>
'..@I=' - ideal fit; '..@V=' - min. significant fitness difference (10^-x)
'..@S=' - max.#stable cycles to stop; '..@E=' - #best models to store <def.50>
'..@P=' - pheromone volatility <def.0.5>; '..@A=' - pheromone impact <def.1>
'..@B=' - model size impact on pheromone update <def.1>
'..@TMIN=,..@TMAX=' - min and max pheromone levels <def. 0.001 and 1>

Pheromone update modes <def. by all ants in proportion to model fitness>:
'..@GLB' - by global best, '..@LOC' - by iter best'
'..@ADG' - extra update by global best at each iter
'..@ASR' - update based on model rank; '..'@PST' - post-optimize



---Particle Swarm settings: (knn+ -O=PS)
'-PS@...' - Particle Swarm settings: e.g. '-PS@N=200@D=100@S=5@V=-1'
'..@N=' - swarm size <def.100>; '..@D=' - max.#cycles <def.300>
'..@S=' - max.#stable cycles to stop; '..@V=' - convergence thershold (10^-x)
'..@W=' - inertia <def.0.7>; '..@VMAX=' - max.velocity <def.10>
'..@B=' - #neighbors to socialize, 0 is for entire swarm <def.0>
'..@C1=,..@C2=' - cognitive and social learning rates  <def. 1.5 and 1.5>

