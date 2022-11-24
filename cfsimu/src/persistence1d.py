import numpy as np


class UnionFind:
    """
    Implements the Union-Find data structure.
    
    It keeps track of elements
    partitioned into a number of disjoint (non-overlapping) sets.
    This data structure is required for merge trees and similar algorithms.
    This implementation uses path compression in several places,
    written with a merge-tree-like algorithm in mind.
    A set is identified by the ID of the defining element (e.g., vertex).
    
    Author: Tino Weinkauf
    """
    NOSET = -1
    
    def __init__(self, NumElements):
        """Initializes the domain with @NumElements elements living in zero sets."""
        self.Domain = [int(UnionFind.NOSET)] * NumElements
        self.NumSets = 0


    def _assert(self, condition, text):
        if not condition:
            raise ValueError(text)
        
        
    def GetNumSets(self):
        """Returns the number of sets."""
        return self.NumSets
    
    
    def MakeSet(self, idxElement):
        """Creates a new set with the given @idxElement as a root."""
        self._assert(self.Domain[idxElement] == UnionFind.NOSET, "You cannot call MakeSet() on an existing set.")
            
        self.Domain[idxElement] = idxElement
        self.NumSets += 1


    def ExtendSet(self, idxElementFrom, idxElementTo):
        """
        Extends a set from one element to the next.
        @note The element identified by @idxElementFrom needs to belong to a set.
        @note The element identified by @idxElementTo cannot belong to a set.
        """
        self._assert(self.Domain[idxElementTo] == UnionFind.NOSET, "You cannot call ExtendSet() on an existing set.")
        self._assert(self.Domain[idxElementFrom] != UnionFind.NOSET, "You cannot call ExtendSet() without an existing set.")
        
        idxRoot = self.Find(idxElementFrom)
        self.Domain[idxElementTo] = idxRoot
        self.Domain[idxElementFrom] = idxRoot #compression only; not functionally required.


    def ExtendSetByID(self, idxRoot, idxElementTo):
        """
        Extends a set with a given set ID.
        @note The set identified by @idxRoot needs to exist.
            In particular, this needs to be true: Find(idxRoot) == idxRoot
        @note The element identified by @idxElementTo cannot belong to a set.
        """
        self._assert(self.Domain[idxElementTo] == UnionFind.NOSET, "You cannot call ExtendSetByID() on an existing set.")
        self._assert(self.Domain[idxRoot] != UnionFind.NOSET, "You cannot call ExtendSetByID() without an existing set.")
        self._assert(idxRoot == self.Find(idxRoot), "Please call the function ExtendSetByID() with an existing set ID.")

        self.Domain[idxElementTo] = idxRoot       


    def Union(self, idxElementMergeThisOne, idxElementIntoThisOne):
        """
        Merges two sets into one.
        The two sets are identified by their elements @idxElementMergeThisOne and @idxElementIntoThisOne.
        The former set is merged into the latter, i.e., the latter one remains.
        This function uses a lot of compression to speed-up later calls to any Find*() function.
        """
        self._assert(idxElementMergeThisOne != UnionFind.NOSET and idxElementIntoThisOne != UnionFind.NOSET, "You can only call Union() on existing sets.")

        #~ Find the root that will survive this union
        SetIDB = self.FindAndCompress(idxElementIntoThisOne)

        #~ Find the root that will merge into the surviving one, compress/merge along the way
        bIsSameSet = self.FindMergeCompressTo(idxElementMergeThisOne, SetIDB)

        #~ Same set?
        if not bIsSameSet: self.NumSets -= 1

        self._assert(self.NumSets > 0, "We should have at least one set.");


    def Find(self, idxElement):
        """
        Finds the ID of the set to which the element @idxElement belongs.
        This function does not use compression, and therefore does not change any underlying data.
        """
        if (self.Domain[idxElement] == UnionFind.NOSET): return UnionFind.NOSET

        #~ Find the root of the set
        while True:
            idxElement = self.Domain[idxElement]
            if (idxElement == self.Domain[idxElement]): break;

        #~ Return root of set == ID of set
        return idxElement


    def FindAndCompress(self, idxElement):
        """
        Finds the ID of the set to which the element @idxElement belongs, and compresses the entire path.
        Compression means that all elements along the path point to the root of the set.
        This makes future calls to any Find*() function faster.
        """
        if (self.Domain[idxElement] == UnionFind.NOSET): return UnionFind.NOSET

        #~ Record a path
        Path = [idxElement];

        #~ Find the root of the set
        while True:
            idxElement = self.Domain[idxElement]
            Path.append(idxElement)
            if (idxElement == self.Domain[idxElement]): break

        #~ Compress the path
        for idxStep in Path: self.Domain[idxStep] = idxElement

        #~ Return root of set == ID of set
        return idxElement


    def FindMergeCompressTo(self, idxElement, idxRoot):
        """
        Find a path from @idxElement to its root and compresses the entire path to a new root.
        Useful only when merging sets.
        @returns true, if the root of @idxElement is already idxRoot, i.e., they belong to the same set.
        @returns false, otherwise.
        """
        if (self.Domain[idxElement] == UnionFind.NOSET): return false

        #~ Record a path
        Path = [idxElement];

        #~ Find the root of the set
        while True:
            idxElement = self.Domain[idxElement]
            Path.append(idxElement)
            if (idxElement == self.Domain[idxElement]): break

        bIsSameSet = (idxElement == idxRoot)

        #~ Compress the path
        for idxStep in Path: self.Domain[idxStep] = idxRoot

        return bIsSameSet
    

def RunPersistence(InputData):
    """
    Finds extrema and their persistence in one-dimensional data.
    
    Local minima and local maxima are extracted, paired,
    and returned together with their persistence.
    The global minimum is extracted as well.
    We assume a connected one-dimensional domain.
    Think of "data on a line", or a function f(x) over some domain xmin <= x <= xmax.
    We are only concerned with the data values f(x)
    and do not care to know the x positions of these values,
    since this would not change which point is a minimum or maximum.
    
    This function returns a list of extrema together with their persistence.
    The list is NOT sorted, but the paired extrema can be identified, i.e.,
    which minimum and maximum were removed together at a particular
    persistence level. As follows:
    The odd entries are minima, the even entries are maxima.
    The minimum at 2*i is paired with the maximum at 2*i+1.
    The last entry of the list is the global minimum.
    It is not paired with a maximum.
    Hence, the list has an odd number of entries.
    
    Author: Tino Weinkauf
    """
 
    #~ How many items do we have?
    NumElements = len(InputData)
    
    #~ Sort data in a stable manner to break ties (leftmost index comes first)
    SortedIdx = np.argsort(InputData, kind='stable')

    #~ Get a union find data structure
    UF = UnionFind(NumElements)

    #~ Paired extrema
    ExtremaAndPersistence = []

    #~ Watershed
    for idx in SortedIdx:
        
        #~ Get neighborhood indices
        LeftIdx = max(idx - 1, 0)
        RightIdx = min(idx + 1, NumElements - 1)
        
        #~ Count number of components in neighborhhood
        NeighborComponents = []
        LeftNeighborComponent = UF.Find(LeftIdx)
        RightNeighborComponent = UF.Find(RightIdx)
        if (LeftNeighborComponent != UnionFind.NOSET): NeighborComponents.append(LeftNeighborComponent)
        if (RightNeighborComponent != UnionFind.NOSET): NeighborComponents.append(RightNeighborComponent)
        #~ Left and Right cannot be the same set in a 1D domain
        #~ self._assert(LeftNeighborComponent == UnionFind.NOSET or RightNeighborComponent == UnionFind.NOSET or LeftNeighborComponent != RightNeighborComponent, "Left and Right cannot be the same set in a 1D domain.")
        NumNeighborComponents = len(NeighborComponents)
        
        if (NumNeighborComponents == 0):
            #~ Create a new component
            UF.MakeSet(idx)
        elif (NumNeighborComponents == 1):
            #~ Extend the one and only component in the neighborhood
            #~ Note that NeighborComponents[0] holds the root of a component, since we called Find() earlier to retrieve it
            UF.ExtendSetByID(NeighborComponents[0], idx)
        else:
            #~ Merge the two components on either side of the current point
            #~ The current point is a maximum. We look for the largest minimum on either side to pair with. That is the smallest hub.
            #~ We look for the lowest minimum first (the one that survives) to break the tie in case of equality: np.argmin returns the first occurence in this case.
            idxLowestNeighborComp = np.argmin(InputData[NeighborComponents])
            idxLowestMinimum = NeighborComponents[idxLowestNeighborComp]
            idxHighestMinimum = NeighborComponents[(idxLowestNeighborComp + 1) % 2]
            UF.ExtendSetByID(idxLowestMinimum, idx)
            UF.Union(idxHighestMinimum, idxLowestMinimum)
            
            #~ Record the two paired extrema: index of minimu, index of maximum, persistence value
            Persistence = InputData[idx] - InputData[idxHighestMinimum]
            ExtremaAndPersistence.append((idxHighestMinimum, Persistence))
            ExtremaAndPersistence.append((idx, Persistence))

    idxGlobalMinimum = UF.Find(0)
    ExtremaAndPersistence.append((idxGlobalMinimum, np.inf))
    #~ print("UF is left with %d sets." % UF.NumSets)
    #~ print("Global minimum at %d with value %g" % (idxGlobalMinimum, InputData[idxGlobalMinimum]))

    return ExtremaAndPersistence