namespace OrganicNamer.Core
{
    public class ElementGraph
    {
        private List<Element> atoms; //list of nodes
        public List<Element> Atoms { get { return atoms; } }
        private int[] ends;
        private List<FunctionalGroup> groups;
        public List<FunctionalGroup> Groups { get { return groups; } }
        private Dictionary<string, (string name, int valency)> periodicTable;

        // Private constructor for FromJsonAtoms
        private ElementGraph()
        {
            atoms = new List<Element>();
            groups = new List<FunctionalGroup>();
        }

        // Helper classes for JSON input

        /// <summary>A single atom in the molecule graph, identified by its index in the request's atom list.</summary>
        public class AtomInput
        {
            /// <summary>Element symbol, e.g. "C", "O", "N", "Cl". Case-sensitive, must match the periodic table.</summary>
            public string Element { get; set; } = "";

            /// <summary>Bonds from this atom to other atoms. Hydrogens may be omitted and will be filled in automatically to satisfy valency.</summary>
            public List<BondInput> Bonds { get; set; } = new();
        }

        /// <summary>A bond from the owning atom to another atom in the same request's atom list.</summary>
        public class BondInput
        {
            /// <summary>Zero-based index into the request's atom list of the atom this bond connects to.</summary>
            public int To { get; set; }

            /// <summary>Bond order: 1 for single, 2 for double, 3 for triple.</summary>
            public int Order { get; set; }
        }

        // Fills in missing hydrogen atoms based on each element's valency deficit.
        // For every atom where sum(bond orders) < valency, the required H atoms are appended.
        // Safe to call on molecules that already have explicit H atoms — those atoms
        // already have full valency so no extra H is added.
        public static List<AtomInput> FillImplicitHydrogens(
            List<AtomInput> atoms,
            Dictionary<string, (string name, int valency)> periodicTable)
        {
            // Clone the list so we never mutate the caller's data
            var result = atoms.Select(a => new AtomInput
            {
                Element = a.Element,
                Bonds = a.Bonds.Select(b => new BondInput { To = b.To, Order = b.Order }).ToList()
            }).ToList();

            int nextIndex = result.Count;
            var newHydrogens = new List<AtomInput>();

            for (int i = 0; i < result.Count; i++)
            {
                if (!periodicTable.ContainsKey(result[i].Element)) continue;
                int valency = periodicTable[result[i].Element].valency;
                int usedValency = result[i].Bonds.Sum(b => b.Order);
                int missingH = valency - usedValency;

                for (int h = 0; h < missingH; h++)
                {
                    result[i].Bonds.Add(new BondInput { To = nextIndex, Order = 1 });
                    newHydrogens.Add(new AtomInput
                    {
                        Element = "H",
                        Bonds = new List<BondInput> { new BondInput { To = i, Order = 1 } }
                    });
                    nextIndex++;
                }
            }

            result.AddRange(newHydrogens);
            return result;
        }

        // Static method to build ElementGraph from JSON atoms
        public static ElementGraph FromJsonAtoms(
            List<AtomInput> atoms,
            Dictionary<string, (string name, int valency)> periodicTable)
        {
            var graph = new ElementGraph();
            graph.periodicTable = periodicTable;

            // Create Element objects from JSON
            foreach (var atom in atoms)
            {
                var (name, valency) = periodicTable[atom.Element];
                // Wave 3 (G2): bond storage must hold ALL declared bonds, not just
                // `valency` of them — a hypervalent nitro N declares 5 bond slots against
                // valency 3 and would overflow AddHalfBond. Over-valence is then rejected
                // deliberately in ValidateValences(), where pattern-consumed atoms are
                // the only exemption.
                int declaredValency = atom.Bonds.Sum(b => b.Order);
                graph.atoms.Add(new Element(atom.Element, name, Math.Max(valency, declaredValency)));
            }

            // Establish bonds from JSON (uses existing AddFullBond method)
            for (int i = 0; i < atoms.Count; i++)
            {
                foreach (var bond in atoms[i].Bonds)
                {
                    if (bond.To > i) // Only process each bond once
                    {
                        graph.AddFullBond(i, bond.To, bond.Order);
                    }
                }
            }

            // Run existing initialization logic
            graph.ends = graph.FindEnds();
            graph.groups = graph.FindGroups();
            graph.RemoveDuplicateGroups();
            graph.ValidateValences(); // Wave 3: must run AFTER FindGroups (needs the PolyatomicGroups)

            return graph;
        }

        public bool IsAnEnd(int index)
        {
            return ends.Contains(index);
        }
        private void AddFullBond(int indexOne, int indexTwo, int bondOrder)
        {
            atoms[indexOne].AddHalfBond(indexTwo, bondOrder);
            atoms[indexTwo].AddHalfBond(indexOne, bondOrder);
        }
        private int PreviousSimilarAtoms(int atomIndex) //how many atoms of same element have already been added
        {
            int count = 1;
            for (int i = 0; i < atomIndex; i++)
            {
                if (atoms[i].Name == atoms[atomIndex].Name)
                {
                    count++;
                }
            }
            return count;
        }
        public int AlkylCounter(int index)
        {
            return atoms[index].BondIndexes.Distinct().Count(x => atoms[x].Name == "Carbon"); //hashset as double bonds are counted once
        }
        public int[] AdjacentAtoms(int index, string ignoreSymbol) //adjacent atoms, H is to be ignored as it is dead end
        {
            return atoms[index].BondIndexes.Distinct().Where(x => atoms[x].Symbol != ignoreSymbol && atoms[x].Symbol != "H").ToArray();
        }
        public int[] AdjacentAtoms(int index)
        {
            return AdjacentAtoms(index, "");
        }
        public int BondOrder(int indexOne, int indexTwo) //not distinct as it returns bond order/weighting
        {
            return atoms[indexOne].BondIndexes.Count(x => x == indexTwo);
        }
        private int[] FindEnds()
        {
            List<int> ends = new List<int>();
            for (int chainIndex = 0; chainIndex < atoms.Count; chainIndex++)
            {
                if (atoms[chainIndex].Name == "Carbon" && !ends.Contains(chainIndex))
                {
                    int alkylCount = AlkylCounter(chainIndex);
                    if (alkylCount == 0 || alkylCount == 1) //isolated carbon ie O-C-O, or unique end
                    {
                        ends.Add(chainIndex);
                    }
                }
            }
            if (ends.Count == 0)
            {
                if (IsCyclic())
                    return Array.Empty<int>(); // signal to IUPAC: this is a ring, not an error
                throw new Exception("The molecule seems to have no ends");
            }
            return ends.ToArray();
        }
        public bool IsCyclic()
        {
            HashSet<int> visited = new HashSet<int>();
            for (int i = 0; i < atoms.Count; i++)
            {
                if (atoms[i].Name == "Carbon" && !visited.Contains(i))
                {
                    if (HasCycleDFS(i, -1, visited))
                        return true;
                }
            }
            return false;
        }
        private bool HasCycleDFS(int current, int parent, HashSet<int> visited)
        {
            visited.Add(current);
            foreach (int neighbour in AdjacentAtoms(current).Where(n => atoms[n].Name == "Carbon"))
            {
                if (!visited.Contains(neighbour))
                {
                    if (HasCycleDFS(neighbour, current, visited))
                        return true;
                }
                else if (neighbour != parent)
                {
                    return true; // back-edge found → cycle
                }
            }
            return false;
        }
        public List<int> FindRing()
        {
            for (int i = 0; i < atoms.Count; i++)
            {
                if (atoms[i].Name != "Carbon" || AlkylCounter(i) < 2)
                    continue;
                int[] cn = AdjacentAtoms(i).Where(n => atoms[n].Name == "Carbon").ToArray();
                for (int a = 0; a < cn.Length; a++)
                {
                    for (int b = a + 1; b < cn.Length; b++)
                    {
                        List<int> partial = FindPathBlocking(cn[a], cn[b], i);
                        if (partial.Count > 0 && partial[partial.Count - 1] == cn[b]) // path really reached cn[b]
                        {
                            partial.Insert(0, i);
                            return partial;
                        }
                    }
                }
            }
            throw new Exception("Cycle detected but no ring could be extracted");
        }
        public List<int> FindPathBlocking(int start, int end, int blocked)
        {
            List<int> path = new List<int>();
            HashSet<int> visited = new HashSet<int> { blocked }; // pre-seed blocked node
            FindPathRecursive(start, end, ref visited, ref path);
            return path;
        }
        public List<(int ringPosition, int atomIndex, bool isCarbon)> GetRingSubstituents(List<int> ring)
        {
            var result = new List<(int, int, bool)>();
            HashSet<int> ringSet = new HashSet<int>(ring);
            for (int pos = 0; pos < ring.Count; pos++)
            {
                int ringC = ring[pos];
                foreach (int neighbour in AdjacentAtoms(ringC))
                {
                    if (!ringSet.Contains(neighbour))
                    {
                        bool isCarbon = atoms[neighbour].Name == "Carbon";
                        result.Add((pos, neighbour, isCarbon));
                    }
                }
            }
            return result;
        }
        public bool IsAromatic(List<int> ring)
        {
            if (ring.Count != 6) return false;

            // All ring atoms must be carbon
            if (ring.Any(i => atoms[i].Name != "Carbon")) return false;

            // Check for alternating bond orders around the ring
            // Try both starting patterns (1,2,1,2,1,2 and 2,1,2,1,2,1)
            for (int startOrder = 1; startOrder <= 2; startOrder++)
            {
                bool matches = true;
                for (int i = 0; i < 6; i++)
                {
                    int expectedOrder = (i % 2 == 0) ? startOrder : (3 - startOrder);
                    int actualOrder = BondOrder(ring[i], ring[(i + 1) % 6]);
                    if (actualOrder != expectedOrder)
                    {
                        matches = false;
                        break;
                    }
                }
                if (matches) return true;
            }
            return false;
        }
        public List<int> CollectReachable(int start, HashSet<int> blocked)
        {
            List<int> result = new List<int>();
            HashSet<int> visited = new HashSet<int>(blocked);
            CollectReachableDFS(start, visited, result);
            return result; // result[0] == start — guaranteed by DFS order; callers rely on this
        }
        private void CollectReachableDFS(int current, HashSet<int> visited, List<int> result)
        {
            visited.Add(current);
            result.Add(current);
            foreach (int neighbour in AdjacentAtoms(current))
            {
                if (!visited.Contains(neighbour))
                {
                    CollectReachableDFS(neighbour, visited, result);
                }
            }
        }
        public List<int> GetBridgingAtoms()
        {
            // FindGroups() already detects chain-breaking heteroatoms: for a
            // CarbonOtherCarbonGroup, MainIndex IS the central bridging atom.
            List<int> bridges = new List<int>();
            foreach (FunctionalGroup g in groups)
            {
                if (g is CarbonOtherCarbonGroup)
                {
                    bridges.Add(g.MainIndex);
                }
            }
            return bridges;
        }
        public (List<int> sideA, List<int> sideB) SplitAtBridgingAtom(int bridgeIndex)
        {
            // AdjacentAtoms already returns distinct, non-H neighbours — reuse it rather
            // than re-filtering BondIndexes (which is -1-padded until H-filling runs).
            int[] carbonNeighbours = AdjacentAtoms(bridgeIndex)
                .Where(n => atoms[n].Name == "Carbon").ToArray();

            if (carbonNeighbours.Length != 2)
                throw new Exception("Bridging atoms must bridge exactly 2 carbons"); // rejects trimethylamine etc.

            var blocked = new HashSet<int> { bridgeIndex };
            List<int> sideA = CollectReachable(carbonNeighbours[0], blocked);

            // Ring-through-heteroatom guard: if the second attachment carbon is reachable
            // from the first WITHOUT the bridge atom, the bridge sits inside a ring
            // (THF, epoxides, oxetane, pyrrolidine, lactones). IsCyclic() cannot see these
            // rings — it only checks the carbon subgraph — so this is the ONLY place they
            // are caught. Without it, both sides come back as the same atom set and the
            // molecule is confidently misnamed (THF → "butoxybutane").
            if (sideA.Contains(carbonNeighbours[1]))
                throw new Exception("Rings containing a heteroatom are not supported");

            List<int> sideB = CollectReachable(carbonNeighbours[1], blocked);
            return (sideA, sideB);
        }
        // True when removing the bridge atom leaves its two attachment carbons still
        // mutually reachable — i.e. the bridge closes a ring rather than joining two
        // separate fragments. Same test as the guard inside SplitAtBridgingAtom, hoisted
        // so the routing decision in IUPAC can consult it before choosing a lane.
        public bool BridgeClosesARing(int bridgeIndex)
        {
            int[] carbonNeighbours = AdjacentAtoms(bridgeIndex)
                .Where(n => atoms[n].Name == "Carbon").ToArray();
            if (carbonNeighbours.Length != 2)
                return false;                     // other guards own the non-2-carbon shapes
            var blocked = new HashSet<int> { bridgeIndex };
            return CollectReachable(carbonNeighbours[0], blocked).Contains(carbonNeighbours[1]);
        }
        // Sibling of SplitAtBridgingAtom for tertiary N-bridges (3 carbon neighbours).
        // The 2-way original stays untouched. Same ring-through-heteroatom guard.
        public List<List<int>> SplitAtBridgingAtomMultiway(int bridgeIndex)
        {
            int[] carbonNeighbours = AdjacentAtoms(bridgeIndex)
                .Where(n => atoms[n].Name == "Carbon").ToArray();

            var blocked = new HashSet<int> { bridgeIndex };
            var sides = new List<List<int>>();
            foreach (int cn in carbonNeighbours)
            {
                List<int> side = CollectReachable(cn, blocked); // side[0] == cn
                foreach (int other in carbonNeighbours)
                    if (other != cn && side.Contains(other))
                        throw new Exception("Rings containing a heteroatom are not supported");
                sides.Add(side);
            }
            return sides;
        }
        public List<List<int>> FindBranches(List<int> path)
        {
            List<List<int>> branches = new List<List<int>>();
            foreach (int i in path)
            {
                if (AlkylCounter(i) > 2)
                {
                    foreach (int end in ends)
                    {
                        if (!path.Contains(end))
                        {
                            List<int> branch = FindPath(i, end);
                            if (branch.Intersect(path).Count() == 1) //new branch
                            {
                                branches.Add(branch);
                            }
                        }
                    }
                }
            }
            return branches;
        }
        public List<List<int>> FindEveryLongestPath()
        {
            if (ends.Length == 1) //only one carbon atom
            {
                return new List<List<int>> { new List<int> { ends[0] } };
            }
            List<List<int>> paths = new List<List<int>>();
            for (int i = 0; i < ends.Length; i++)
            {
                for (int j = 0; j < ends.Length; j++)
                {
                    if (i != j)
                    {
                        paths.Add(FindPath(ends[i], ends[j]));
                    }

                }
            }
            return paths;
        }
        public List<int> FindPath(int start, int end)
        {
            List<int> path = new List<int>();
            HashSet<int> visited = new HashSet<int>();
            FindPathRecursive(start, end, ref visited, ref path); //ref to keep track through the whole stack
            return path;
        }
        private bool FindPathRecursive(int current, int end, ref HashSet<int> visited, ref List<int> path) //recurisve DFS
        {
            visited.Add(current);
            path.Add(current);
            if (current == end)
            {
                return true;
            }
            int[] bonds = AdjacentAtoms(current);
            foreach (int bond in bonds)
            {
                if (!visited.Contains(bond))
                {
                    if (FindPathRecursive(bond, end, ref visited, ref path)) //found end somewhere down stack
                    {
                        return true;
                    }
                }
            }
            path.Remove(current); //hit deadend somewhere down stack. backtrack
            return false;
        }
        // G2 (Wave 3): pattern-based polyatomic group recognition, run BEFORE the
        // generic per-atom loop in FindGroups. One pattern so far: nitro, accepted
        // ONLY as the hypervalent neutral form N(=O)(=O) - one N with exactly one
        // single-bonded C, exactly two double-bonded TERMINAL O, and nothing else.
        // Matched heteroatoms are returned as the consumed set; the generic loop
        // (including the heteroatom-heteroatom throw) skips them. Every non-matching
        // shape (charge-separated nitro, nitroso, nitrite/nitrate esters, N-N,
        // peroxides) is not consumed and still throws exactly as before.
        private HashSet<int> FindPolyatomicGroups(List<FunctionalGroup> newGroups)
        {
            HashSet<int> consumed = new HashSet<int>();
            for (int n = 0; n < atoms.Count; n++)
            {
                if (atoms[n].Symbol != "N")
                    continue;
                int[] neighbours = atoms[n].BondIndexes.Where(x => x >= 0).Distinct().ToArray();
                if (neighbours.Length != 3) //must be exactly C + O + O (an H neighbour kills the match)
                    continue;
                int[] carbons = neighbours.Where(x => atoms[x].Name == "Carbon").ToArray();
                int[] oxygens = neighbours.Where(x => atoms[x].Name == "Oxygen").ToArray();
                if (carbons.Length != 1 || oxygens.Length != 2)
                    continue;
                if (BondOrder(n, carbons[0]) != 1)
                    continue;
                if (oxygens.Any(o => BondOrder(n, o) != 2))
                    continue;
                // Terminal-O check: each O bonds to nothing but this N. Not implied by
                // valence (input valence is unvalidated here) - an adversarial O
                // bridging onward must not be consumed.
                if (oxygens.Any(o => atoms[o].BondIndexes.Where(x => x >= 0).Any(x => x != n)))
                    continue;
                newGroups.Add(new PolyatomicGroup("NO2", carbons[0], new[] { n, oxygens[0], oxygens[1] }));
                consumed.Add(n);
                consumed.Add(oxygens[0]);
                consumed.Add(oxygens[1]);
            }
            return consumed;
        }
        private List<FunctionalGroup> FindGroups()
        {
            List<FunctionalGroup> newGroups = new List<FunctionalGroup>();
            HashSet<int> consumed = FindPolyatomicGroups(newGroups); //G2 pattern pass (Wave 3)
            for (int atomIndex = 0; atomIndex < atoms.Count; atomIndex++)
            {
                Element atom = atoms[atomIndex];
                if (atom.Name == "Carbon")
                {
                    int[] unvisitedBonds = AdjacentAtoms(atomIndex);
                    foreach (int bondIndex in unvisitedBonds)
                    {
                        if (consumed.Contains(bondIndex))
                            continue; //pattern atom (nitro N): suppress the loose C-N group (R4)
                        Element bondedAtom = atoms[bondIndex];
                        int order = BondOrder(atomIndex, bondIndex);
                        if (bondedAtom.Name == "Carbon" && order != 1) //two carbons forming double / triple bond
                        {
                            newGroups.Add(new CarbonCarbonGroup(order, atomIndex, bondIndex)); //
                        }
                        else if (bondedAtom.Name != "Carbon") //normal group eg C=O, C-N
                        {
                            newGroups.Add(new FunctionalGroup(bondedAtom.Symbol, order, atomIndex));
                        }
                    }
                }
                else if (consumed.Contains(atomIndex))
                {
                    //member of a matched polyatomic pattern: fully described by its
                    //PolyatomicGroup - skip the heteroatom-heteroatom throw and the
                    //bridging-atom detection (R4)
                }
                else if (atom.Name != "Hydrogen" && AdjacentAtoms(atomIndex, "C").Count() > 0)  //non-carbon atom bonded to some non-carbon atoms
                {
                    throw new Exception("Multiple non-carbon atoms are bonded together");
                }
                else if (AlkylCounter(atomIndex) > 1) //C-O-C
                {
                    List<int> carbonIndexes = new List<int>();
                    foreach (int bond in atom.BondIndexes)
                    {
                        if (atoms[bond].Name == "Carbon")
                        {
                            if (BondOrder(atomIndex, bond) > 1)
                            {
                                throw new Exception("An atom makes non-single bonds with multiple carbons"); //cannot have CarbonOthercarbon groups where the atom makes double/triple bonds with carbon
                            }
                            carbonIndexes.Add(bond);
                        }
                    }
                    newGroups.Add(new CarbonOtherCarbonGroup(atom.Symbol, atomIndex, carbonIndexes.ToArray()));
                }
            }
            return newGroups;
        }
        // Wave 3 (G2): bond storage is no longer capped at valency (see FromJsonAtoms),
        // so over-valent input must be rejected explicitly instead of crashing the
        // array. Atoms consumed as members of a matched polyatomic pattern (the
        // hypervalent nitro N) are exempt - the pattern already pinned their bonding.
        private void ValidateValences()
        {
            HashSet<int> patternMembers = new HashSet<int>(
                groups.OfType<PolyatomicGroup>().SelectMany(g => g.MemberIndexes));
            for (int i = 0; i < atoms.Count; i++)
            {
                int usedBonds = atoms[i].BondIndexes.Count(x => x >= 0);
                if (usedBonds > periodicTable[atoms[i].Symbol].valency && !patternMembers.Contains(i))
                    throw new Exception($"Atom {i} ({atoms[i].Symbol}) has more bonds than its valency allows");
            }
        }
        private void RemoveDuplicateGroups() //CarbonCarbonGroups are recorded twice (once for each carbon)
        {
            for (int i = groups.Count - 1; i >= 0; i--)
            {
                if (groups[i] is CarbonCarbonGroup)
                {
                    CarbonCarbonGroup cgroupOne = (CarbonCarbonGroup)groups[i];
                    for (int j = groups.Count - 1; j >= 0; j--)
                    {
                        if (groups[j] is CarbonCarbonGroup && i != j)
                        {
                            CarbonCarbonGroup cgroupTwo = (CarbonCarbonGroup)groups[j];
                            if (cgroupOne.IsSame(cgroupTwo))
                            {
                                groups.RemoveAt(i);
                            }
                        }
                    }
                }
            }
        }
        public void MergeFunctionalGroups(Dictionary<HashSet<string>, string> toMerge)
        {
            for (int i = 0; i < atoms.Count; i++)
            {
                HashSet<string> formulae = groups.Where(x => x.Involves(i)).Select(x => x.GroupFormula).ToHashSet();
                foreach (HashSet<string> collect in toMerge.Keys)
                {
                    if (collect.IsSubsetOf(formulae)) //enough groups to merge
                    {
                        formulae = new HashSet<string>(collect); //new set, otherwise it is passsed by ref. and collect is mutated
                        groups.Add(new MergedGroup(toMerge[collect], i));
                        for (int j = groups.Count - 1; j >= 0; j--) //negative iteration as mutating data structure
                        {
                            if (groups[j].Involves(i) && formulae.Contains(groups[j].GroupFormula))
                            {
                                formulae.Remove(groups[j].GroupFormula);
                                groups.RemoveAt(j);
                            }
                        }
                    }
                }
            }
        }
        public static void SavePeriodicTable(string fileName, Dictionary<string, (string name, int valency)> newPeriodicTable)
        {
            if (!fileName.EndsWith(".PeriodicTable"))
            {
                fileName += ".PeriodicTable";
            }
            using (BinaryWriter writefile = new BinaryWriter(File.Open(fileName, FileMode.Create)))
            {
                writefile.Write(Convert.ToInt16(newPeriodicTable.Count));
                foreach (KeyValuePair<string, (string name, int valency)> element in newPeriodicTable)
                {
                    writefile.Write(element.Key); //must contain a "Carbon" "C" and "Hydrogen" "H"
                    writefile.Write(element.Value.name);
                    writefile.Write(Convert.ToInt16(element.Value.valency));
                }
                writefile.Close();
            }
        }
        public static void RestoreDefaultPeriodicTable()
        {
            Dictionary<string, (string name, int valency)> newPeriodicTable = new Dictionary<string, (string name, int valency)> {
                { "C", ("Carbon", 4) },
                { "H", ("Hydrogen", 1) },
                { "O", ("Oxygen", 2) },
                { "N", ("Nitrogen", 3) },
                { "F", ("Fluorine", 1) },
                { "Cl", ("Chlorine", 1) },
                { "Br", ("Bromine", 1) },
                { "I", ("Iodine", 1) },
                { "S", ("Sulphur", 2) } };
            SavePeriodicTable("Default", newPeriodicTable);
        }
        private Dictionary<string, (string name, int valency)> LoadPeriodicTable(string fileName)
        {
            if (!fileName.EndsWith(".PeriodicTable"))
            {
                fileName += ".PeriodicTable";
            }
            Dictionary<string, (string name, int valency)> newPeriodicTable = new Dictionary<string, (string name, int valency)>();
            using (BinaryReader readfile = new BinaryReader(File.Open(fileName, FileMode.Open)))
            {
                int elementCount = readfile.ReadInt16();
                for (int i = 0; i < elementCount; i++)
                {
                    string symbol = readfile.ReadString();
                    string name = readfile.ReadString();
                    int valency = readfile.ReadInt16();
                    newPeriodicTable.Add(symbol, (name, valency));
                }
                readfile.Close();
            }
            return newPeriodicTable;
        }
    }
}