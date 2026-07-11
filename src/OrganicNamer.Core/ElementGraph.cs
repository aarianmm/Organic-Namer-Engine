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
        public class AtomInput
        {
            public string Element { get; set; } = "";
            public List<BondInput> Bonds { get; set; } = new();
        }

        public class BondInput
        {
            public int To { get; set; }
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
                graph.atoms.Add(new Element(atom.Element, name, valency));
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
            // Use a carbon with exactly 2 C-neighbours (pure ring carbon) as anchor
            int ringCarbon = -1;
            for (int i = 0; i < atoms.Count; i++)
            {
                if (atoms[i].Name == "Carbon" && AlkylCounter(i) == 2)
                {
                    ringCarbon = i;
                    break;
                }
            }
            if (ringCarbon == -1) // fallback: all ring carbons have substituents
            {
                for (int i = 0; i < atoms.Count; i++)
                {
                    if (atoms[i].Name == "Carbon" && AlkylCounter(i) >= 2)
                    {
                        ringCarbon = i;
                        break;
                    }
                }
            }
            int[] carbonNeighbours = AdjacentAtoms(ringCarbon)
                .Where(n => atoms[n].Name == "Carbon").ToArray();
            int c1 = carbonNeighbours[0];
            int c2 = carbonNeighbours[1];
            List<int> partialRing = FindPathBlocking(c1, c2, ringCarbon);
            partialRing.Insert(0, ringCarbon);
            return partialRing;
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
        private List<FunctionalGroup> FindGroups()
        {
            List<FunctionalGroup> newGroups = new List<FunctionalGroup>();
            for (int atomIndex = 0; atomIndex < atoms.Count; atomIndex++)
            {
                Element atom = atoms[atomIndex];
                if (atom.Name == "Carbon")
                {
                    int[] unvisitedBonds = AdjacentAtoms(atomIndex);
                    foreach (int bondIndex in unvisitedBonds)
                    {
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