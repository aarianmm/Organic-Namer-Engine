using System.Text.RegularExpressions;

namespace OrganicNamer.Core
{
    public class IUPAC
    {
        private readonly List<FunctionalGroup> groups;
        private readonly ElementGraph atoms;
        public struct namingSpec
        {
            public string[] alkylNames;
            public string[] numericalPrefixes;
            public Dictionary<HashSet<string>, string> merging; //some groups are combined and named as one ie C-O and C=O to COOH
            public Dictionary<string, string> prefixOnly;
            public Dictionary<string, (string name, int priority)> middle;
            public Dictionary<string, (string prefix, string suffix, int priority)> prefixOrSuffix;
            public Dictionary<string, ((string prefix, string suffix, int priority) middle, (string prefix, string suffix, int priority) end)> endDependentPrefixOrSuffix;
        }
        private readonly namingSpec spec;
        public string[] names;
        private List<(List<int> chain, List<List<int>> branches)> allChainsAndBranches;
        private string suffixFormula;
        private string suffixRoot;
        private bool suffixIsEnd;
        private bool suffixIsMiddle; //if the position doesnt matter, both booleans are false 

        public IUPAC(namingSpec spec, ElementGraph atoms)
        {
            this.spec = spec;
            this.atoms = atoms;
            groups = atoms.Groups;

            // ──── Guards (Phase 3) ────
            // MUST be before MergeFunctionalGroups — the merge step would
            // incorrectly combine C=O + C-O into COOH on the ester carbon.
            List<int> bridgingAtoms = atoms.GetBridgingAtoms();
            bool cyclic = atoms.IsCyclic();

            // A molecule satisfying both early-exit conditions (anisole,
            // methoxycyclohexane) would be confidently misnamed by whichever
            // exit ran first — neither path can name it, so reject up front.
            if (bridgingAtoms.Count > 0 && cyclic)
                throw new Exception("Molecules with both a ring and a bridging heteroatom are not supported");
            if (bridgingAtoms.Count > 1)
                throw new Exception("Multiple bridging heteroatoms are not supported");

            // Early exit: bridging molecules (esters, ethers, amines)
            if (bridgingAtoms.Count == 1)
            {
                names = NameBridgedMolecule(bridgingAtoms[0]);
                names = names.Distinct().ToArray();
                suffixFormula = "";
                suffixRoot = "";
                return;
            }

            // Early exit: cyclic molecules (rings) — must be before FindEveryLongestPath
            // which requires a non-empty ends array.
            if (cyclic)
            {
                names = NameCyclicMolecule();
                names = names.Distinct().ToArray();
                suffixFormula = "";
                suffixRoot = "";
                return;
            }

            // Standard chain-based path (unchanged logic)
            atoms.MergeFunctionalGroups(spec.merging);
            CheckGroups(groups); //make sure the binary file conatins info for all the groups in this molecule
            //naming
            List<List<int>> possibleChains = atoms.FindEveryLongestPath();
            allChainsAndBranches = NarrowDownChainsByBranches(possibleChains);
            suffixFormula = "";
            suffixRoot = "";
            FindHighestPrioritySuffix();
            if (possibleChains.Count == 0)
            {
                throw new Exception("Impossibe to name, as branches must be empty.");
            }
            if (possibleChains.Count != 1)
            {
                NarrowDownChainsByLength();
                if (possibleChains.Count != 1)
                {
                    NarrowDownChainsBySuffix();
                    if (possibleChains.Count != 1)
                    {
                        NarrowDownChainsByMiddle();
                        if (possibleChains.Count != 1)
                        {
                            NarrowDownChainsByPrefixes();
                        }
                    }
                }
            }
            names = ConstructName();
            names = names.Distinct().ToArray();
        }
        private string[] ConstructName()
        {
            string[] names = new string[allChainsAndBranches.Count];
            for (int i = 0; i < allChainsAndBranches.Count; i++)
            {
                List<int> chain = allChainsAndBranches[i].chain;
                List<List<int>> branches = allChainsAndBranches[i].branches;
                string longestAlkylName = spec.alkylNames[chain.Count];
                Dictionary<string, List<int>> prefixCarbons = FindPrefixCarbons(chain, branches);
                Dictionary<string, List<int>> middleCarbons = FindMiddleCarbons(chain);
                List<int> suffixCarbonsList = FindSuffixCarbons(chain);
                // If no suffix groups are on this chain, suppress the suffix (the group is on a branch)
                Dictionary<string, List<int>> suffixCarbons = (suffixCarbonsList.Count == 0 && suffixRoot != "")
                    ? new Dictionary<string, List<int>> { { "", new List<int>() } }
                    : new Dictionary<string, List<int>> { { suffixRoot, suffixCarbonsList } };
                string prefixName = NameSegment(prefixCarbons);
                string middleName = NameSegment(middleCarbons);
                string suffixName = NameSegment(suffixCarbons);
                names[i] = FormatName(prefixName + longestAlkylName + middleName + suffixName);
            }
            return names;
        }
        private string FormatName(string name)
        {
            name = name.ToLower();
            Regex vowels = new Regex(@"\|(.)(?=[^a-z]*[aeiouy])"); //optional vowel followed by vowel is removed
            Regex startDashes = new Regex(@"([a-z])(\d)");
            Regex endDashes = new Regex(@"(\d)([a-z])");
            Regex parenthesisDashes = new Regex(@"(\d)(\()"); //digit before parenthesis needs a dash
            Regex closeParenDashes = new Regex(@"(\))(\d)");  //digit after a closing parenthesis too
            name = vowels.Replace(name, "");
            name = name.Replace("|", "");
            name = startDashes.Replace(name, (m) => m.Groups[1].Value + "-" + m.Groups[2].Value);
            name = endDashes.Replace(name, (m) => m.Groups[1].Value + "-" + m.Groups[2].Value);
            name = parenthesisDashes.Replace(name, (m) => m.Groups[1].Value + "-" + m.Groups[2].Value);
            name = closeParenDashes.Replace(name, (m) => m.Groups[1].Value + "-" + m.Groups[2].Value);
            return name;
        }
        private List<(List<int> chain, List<List<int>> branches)> NarrowDownChainsByBranches(List<List<int>> chains)
        {
            List<(List<int> chain, List<List<int>> branches)> allChainsAndBranches = new List<(List<int> chain, List<List<int>> branches)>();
            foreach (List<int> chain in chains)
            {
                List<List<int>> branches = atoms.FindBranches(chain);
                allChainsAndBranches.Add((chain, branches));
            }
            for (int i = allChainsAndBranches.Count - 1; i >= 0; i--)
            {
                List<int> chain = allChainsAndBranches[i].chain;
                List<List<int>> branches = allChainsAndBranches[i].branches;
                if (!CheckBranchValidity(chain, branches))
                {
                    allChainsAndBranches.RemoveAt(i); //branch invalid, so whole chain invalid
                }
            }
            return allChainsAndBranches;
        }
        private void NarrowDownChainsByPrefixes()
        {
            List<int> prefixesCarbonSums = new List<int>();
            for (int i = 0; i < allChainsAndBranches.Count; i++)
            {
                List<int> chain = allChainsAndBranches[i].chain;
                List<List<int>> branches = allChainsAndBranches[i].branches;
                List<List<int>> prefixCarbons = FindPrefixCarbons(chain, branches).Values.ToList();
                int prefixCarbonSum = 0;
                foreach (List<int> a in prefixCarbons)
                {
                    prefixCarbonSum += a.Sum();
                }
                prefixesCarbonSums.Add(prefixCarbonSum);
            }
            int lowestSum = prefixesCarbonSums.Min();
            for (int i = allChainsAndBranches.Count - 1; i >= 0; i--)  //negative iter. as mutating size of list
            {
                if (prefixesCarbonSums[i] != lowestSum)
                {
                    allChainsAndBranches.RemoveAt(i);
                }
            }
        }
        private void NarrowDownChainsByMiddle()
        {
            string middlePriorityFormula = FindHighestPriorityMiddleFormula();
            string middlePriorityName = spec.middle[middlePriorityFormula].name;
            List<int> middleCarbonSums = new List<int>();
            for (int i = 0; i < allChainsAndBranches.Count; i++)
            {
                List<int> chain = allChainsAndBranches[i].chain;
                List<int> middleCarbons = FindMiddleCarbons(chain)[middlePriorityName];
                int middleCarbonSum = middleCarbons.Sum();
                middleCarbonSums.Add(middleCarbonSum);
            }
            int lowestSum = middleCarbonSums.Min();
            for (int i = allChainsAndBranches.Count - 1; i >= 0; i--)
            {
                if (middleCarbonSums[i] != lowestSum)
                {
                    allChainsAndBranches.RemoveAt(i);
                }
            }
        }
        private void NarrowDownChainsBySuffix()
        {
            List<int> suffixCarbonSums = new List<int>();
            for (int i = 0; i < allChainsAndBranches.Count; i++)
            {
                List<int> chain = allChainsAndBranches[i].chain;
                List<int> suffixCarbons = FindSuffixCarbons(chain);
                int suffixCarbonSum = suffixCarbons.Sum();
                suffixCarbonSums.Add(suffixCarbonSum);
            }
            int lowestSum = suffixCarbonSums.Min();
            for (int i = allChainsAndBranches.Count - 1; i >= 0; i--)
            {
                if (suffixCarbonSums[i] != lowestSum)
                {
                    allChainsAndBranches.RemoveAt(i);
                }
            }
        }
        private void NarrowDownChainsByLength()
        {
            int longest = allChainsAndBranches.MaxBy(x => x.chain.Count).chain.Count;
            for (int i = allChainsAndBranches.Count - 1; i >= 0; i--)
            {
                if (allChainsAndBranches[i].chain.Count < longest)
                {
                    allChainsAndBranches.RemoveAt(i);
                }
            }
        }
        private Dictionary<string, List<int>> FindPrefixCarbons(List<int> chain, List<List<int>> branches)
        {
            Dictionary<string, List<int>> prefixesAndIndexes = new Dictionary<string, List<int>>();
            foreach (FunctionalGroup group in groups)
            {
                if (chain.IndexOf(group.MainIndex) == -1)
                    continue; // skip groups on branch atoms — handled by BuildSubstituentName

                string name = "";
                int carbonNumber = chain.IndexOf(group.MainIndex) + 1; //carbon number is position in carbon chain
                if (group.GroupFormula != suffixFormula) //group is not the suffix one, so belongs in the prefix
                {
                    if (spec.prefixOnly.ContainsKey(group.GroupFormula))
                    {
                        name = spec.prefixOnly[group.GroupFormula];
                    }
                    else if (spec.prefixOrSuffix.ContainsKey(group.GroupFormula))
                    {
                        name = spec.prefixOrSuffix[group.GroupFormula].prefix;
                    }
                }
                if (spec.endDependentPrefixOrSuffix.ContainsKey(group.GroupFormula))
                {
                    if (atoms.IsAnEnd(group.MainIndex) && !(suffixIsEnd && group.GroupFormula == suffixFormula)) //group is on the end and isnt the same as the suffix one
                    {
                        name = spec.endDependentPrefixOrSuffix[group.GroupFormula].end.prefix;
                    }
                    else if (!atoms.IsAnEnd(group.MainIndex) && !(suffixIsMiddle && group.GroupFormula == suffixFormula)) //group is in the middle and isnt the same as the suffix one
                    {
                        name = spec.endDependentPrefixOrSuffix[group.GroupFormula].middle.prefix;
                    }
                }
                if (name != "")
                {
                    if (prefixesAndIndexes.ContainsKey(name))
                    {
                        prefixesAndIndexes[name].Add(carbonNumber);
                    }
                    else
                    {
                        prefixesAndIndexes.Add(name, new List<int> { carbonNumber });
                    }
                }
            }
            foreach (List<int> branch in branches)
            {
                int carbonNumber = chain.IndexOf(branch[0]) + 1;
                string name = BuildSubstituentName(chain, branch);
                if (prefixesAndIndexes.ContainsKey(name))
                {
                    prefixesAndIndexes[name].Add(carbonNumber);
                }
                else
                {
                    prefixesAndIndexes.Add(name, new List<int> { carbonNumber });
                }
            }
            return prefixesAndIndexes;
        }
        private Dictionary<string, List<int>> FindMiddleCarbons(List<int> chain)
        {
            Dictionary<string, List<int>> middleAndIndexes = new Dictionary<string, List<int>>();
            foreach (FunctionalGroup group in groups)
            {
                if (group is CarbonCarbonGroup)
                {
                    CarbonCarbonGroup cgroup = (CarbonCarbonGroup)group;
                    int carbonNumber;
                    int carbonNumberOne = chain.IndexOf(cgroup.MainIndex) + 1;
                    int carbonNumberTwo = chain.IndexOf(cgroup.OtherCarbonIndex) + 1;
                    if (carbonNumberOne < carbonNumberTwo) //these groups contain two carbons, so have two numbers to choose from. the lowest is chosen
                    {
                        carbonNumber = carbonNumberOne;
                    }
                    else
                    {
                        carbonNumber = carbonNumberTwo;
                    }
                    string name = spec.middle[group.GroupFormula].name;
                    if (middleAndIndexes.ContainsKey(name))
                    {
                        middleAndIndexes[name].Add(carbonNumber);
                    }
                    else
                    {
                        middleAndIndexes.Add(name, new List<int> { carbonNumber });
                    }
                }
            }
            if (middleAndIndexes.Count == 0)
            {
                middleAndIndexes.Add(spec.middle[""].name, new List<int>());
            }
            return middleAndIndexes;
        }
        private List<int> FindSuffixCarbons(List<int> chain)
        {
            List<int> indexes = new List<int>();
            foreach (FunctionalGroup group in groups)
            {
                if (chain.IndexOf(group.MainIndex) == -1)
                    continue; // skip groups on branch atoms

                bool endsAreInvolved = suffixIsMiddle || suffixIsEnd;
                bool groupAndSuffixAreMiddle = suffixIsMiddle && !atoms.IsAnEnd(group.MainIndex);
                bool groupAndSuffixAreEnd = suffixIsEnd && atoms.IsAnEnd(group.MainIndex);
                if (group.GroupFormula == suffixFormula && (!endsAreInvolved || groupAndSuffixAreMiddle || groupAndSuffixAreEnd))
                {
                    indexes.Add(chain.IndexOf(group.MainIndex) + 1);
                }
            }
            return indexes;
        }
        // IUPAC alphabetisation key: compare complete substituent names by their letters
        // only — "(1-methylethyl)" sorts under "m", "(hydroxymethyl)" under "h". Also
        // neutralises the '|' elision markers in raw spec words ("meth|ayl" → "methayl").
        private static string AlphaKey(string name) =>
            new string(name.Where(char.IsLetter).ToArray());
        private string NameSegment(Dictionary<string, List<int>> namesAndCarbonNumbers)
        {
            List<(string numbers, string name)> names = new List<(string, string)>();
            foreach (KeyValuePair<string, List<int>> groupData in namesAndCarbonNumbers)
            {
                string name = groupData.Key;
                List<int> carbonNumbers = groupData.Value;
                carbonNumbers.Sort();
                string numericalPrefix = spec.numericalPrefixes[carbonNumbers.Count];
                string numbers = "";
                for (int i = 0; i < carbonNumbers.Count; i++)
                {
                    numbers += carbonNumbers[i];
                    if (i != carbonNumbers.Count - 1)
                    {
                        numbers += ",";
                    }
                }
                numbers += numericalPrefix;
                names.Add((numbers, name));
            }
            names = names.OrderBy(x => AlphaKey(x.name), StringComparer.Ordinal)
                         .ThenBy(x => x.name, StringComparer.Ordinal).ToList();
            string nameSegment = "";
            for (int i = 0; i < names.Count; i++)
            {
                nameSegment += names[i].numbers + names[i].name;
            }
            return nameSegment;
        }
        private bool CheckBranchValidity(List<int> chain, List<List<int>> branches)
        {
            foreach (List<int> branch in branches)
            {
                for (int i = 1; i < branch.Count; i++)
                {
                    // Sub-branches off branches: not supported (D9 — FindBranches emits
                    // overlapping paths for these; must stay BEFORE the G1 try below).
                    if (atoms.AlkylCounter(branch[i]) > 2)
                        return false;

                    // Groups on intermediate (non-tip) branch carbons: not supported
                    bool isAtTip = (i == branch.Count - 1);
                    if (!isAtTip && groups.Any(g => g.Involves(branch[i])))
                        return false;
                }
                // Wave 2 (D8): a chain candidate is only valid if every branch is
                // G1-nameable. Filtering here (not throwing later) lets a nameable
                // sibling chain win instead of a garbage prefix name.
                try { BuildSubstituentName(chain, branch); }
                catch { return false; }
            }
            return true;
        }
        // ── Phase 4 → Wave 2: chain-branch adapter over G1 ────────────────────────
        // FindBranches produces PATHS (branch[0] = the junction carbon ON the main
        // chain); G1 wants the substituent's own first atom plus a blocked set.
        private string BuildSubstituentName(List<int> chain, List<int> branch)
        {
            return NameSubstituent(branch[0], branch[1], new HashSet<int>(chain),
                                   allowSuffixCapableTip: true);
        }

        // ── Phase 3 / Wave 1: Bridged-molecule naming ──────────────────────────────
        // (esters, ethers, secondary amines, anhydrides, secondary/tertiary amides)
        private enum SideKind { PlainAlkyl, AcidSide }   // E6 adds: AromaticRing
        private enum AcidClass { Ester, Amide, Anhydride }

        private readonly struct SideInfo
        {
            public SideKind Kind { get; init; }
            public int CarbonCount { get; init; }
            // E6 will add ring-related fields here (e.g. the ring atom list) without
            // touching any caller that only reads Kind / CarbonCount.
        }

        private string[] NameBridgedMolecule(int bridgeIndex)
        {
            string bridgeSymbol = atoms.Atoms[bridgeIndex].Symbol;
            int[] carbonNeighbours = atoms.AdjacentAtoms(bridgeIndex)
                                          .Where(n => atoms.Atoms[n].Name == "Carbon").ToArray();
            List<int> acidCarbons = FindAcidCarbons(bridgeIndex);

            // Tertiary N-bridge (3 carbon neighbours): amide only. Checked BEFORE
            // SplitAtBridgingAtom, which rejects ≠2 neighbours.
            if (bridgeSymbol == "N" && carbonNeighbours.Length == 3)
            {
                if (acidCarbons.Count == 1)
                    return NameTertiaryAmide(bridgeIndex);
                throw new Exception("Tertiary amines and N-centred imides are not supported");
            }

            var (sideA, sideB) = atoms.SplitAtBridgingAtom(bridgeIndex); // exactly-2 + ring guard

            if (bridgeSymbol == "O")
                switch (acidCarbons.Count)
                {
                    case 0: return NameEther(sideA, sideB);
                    case 1: return NameEster(sideA, sideB);
                    case 2: return NameAnhydride(sideA, sideB);
                }
            else if (bridgeSymbol == "N")
                switch (acidCarbons.Count)
                {
                    case 0: return NameSecondaryAmine(sideA, sideB);
                    case 1: return NameSecondaryAmide(sideA, sideB);
                    // 2 → imide: rejected by the fall-through throw.
                }

            throw new Exception(
                $"Bridging atom '{bridgeSymbol}' with {acidCarbons.Count} adjacent carbonyl carbon(s) is not supported");
        }

        // Generalises the old FindAcidCarbon (rename + return all, not just the first).
        private List<int> FindAcidCarbons(int bridgeIndex)
        {
            List<int> acidCarbons = new List<int>();
            foreach (int c in atoms.AdjacentAtoms(bridgeIndex)
                                   .Where(n => atoms.Atoms[n].Name == "Carbon"))
                if (groups.Any(g => g.GroupFormula == "C=O" && g.Involves(c)))
                    acidCarbons.Add(c);
            return acidCarbons;
        }

        // Validate a bridged side and classify it. Self-detects an acid side (its
        // bridge-attachment carbon carries a C=O) so callers no longer pass an acid-carbon
        // index or an allowedOxygen. Throws on anything not yet nameable — the capability
        // gate. E6 inserts an aromatic-ring branch here (returning SideKind.AromaticRing)
        // before the plain-alkyl validation loop.
        private SideInfo ClassifySide(List<int> side)
        {
            int attach = side[0]; // bridge-attachment atom (== the bridge's carbon neighbour)

            bool isAcid = atoms.Atoms[attach].Name == "Carbon"
                          && groups.Any(g => g.GroupFormula == "C=O" && g.Involves(attach));
            int allowedOxygen = isAcid ? FindCarbonylOxygen(attach) : -1;

            // Plain-alkyl validation — byte-identical to the old ValidatePlainAlkylSide loop.
            foreach (int idx in side)
            {
                if (idx == allowedOxygen)
                    continue;
                if (atoms.Atoms[idx].Name != "Carbon")
                    throw new Exception("Bridged molecules with additional functional groups are not supported");
                if (groups.Any(g => g is CarbonCarbonGroup && g.Involves(idx)))
                    throw new Exception("Unsaturated bridged chains are not supported");
                int limit = (idx == attach) ? 1 : 2; // attachment carbon must be a chain end
                if (atoms.AlkylCounter(idx) > limit)
                    throw new Exception("Branched or mid-chain-attached bridged chains are not supported");
            }

            int carbonCount = side.Count(i => atoms.Atoms[i].Name == "Carbon");
            return new SideInfo
            {
                Kind = isAcid ? SideKind.AcidSide : SideKind.PlainAlkyl,
                CarbonCount = carbonCount
            };
        }

        // The single carbonyl O double-bonded to an acid carbon. The bridge O (also an
        // oxygen neighbour) is excluded by BondOrder == 2 — bridge bonds are always single.
        private int FindCarbonylOxygen(int acidCarbon) =>
            atoms.AdjacentAtoms(acidCarbon)
                 .First(n => atoms.Atoms[n].Name == "Oxygen" && atoms.BondOrder(acidCarbon, n) == 2);

        // Name a side used as an alkyl substituent / base ("methyl", "ethyl").
        // E6: case SideKind.AromaticRing => "phenyl" (or the substituted-ring name).
        private string AlkylSideName(SideInfo side) => side.Kind switch
        {
            SideKind.PlainAlkyl => FormatName(spec.alkylNames[side.CarbonCount] + "yl"),
            _ => throw new Exception("This side cannot be named as an alkyl substituent")
        };

        // Name an acid side for a given functional class ("ethanoate"/"ethanamide"/"ethanoic").
        // The class owns its suffix; E6 teaches this method the retained aromatic forms
        // (benzoate / benzamide / benzoic) for an acid side whose R is a ring.
        private string AcidSideName(SideInfo side, AcidClass cls)
        {
            string suffix = cls switch
            {
                AcidClass.Ester => "anoate",       // NOT "an|oate": elision would strip the o. §5
                AcidClass.Amide => "anamide",      // NOT "an|amide"
                AcidClass.Anhydride => "anoic",    // NOT "an|oic"
                _ => throw new Exception("Unknown acid class")
            };
            return side.Kind switch
            {
                SideKind.AcidSide => FormatName(spec.alkylNames[side.CarbonCount] + suffix),
                _ => throw new Exception("This side cannot be named as an acid side")
            };
        }

        // Build the "N-…" locant scaffold. Uppercase N survives because each substituent word
        // is FormatName'd by the caller (AlkylSideName) before it arrives.
        //   ["methyl"]           -> "N-methyl"
        //   ["methyl","methyl"]  -> "N,N-dimethyl"
        //   ["ethyl","methyl"]   -> "N-ethyl-N-methyl"  (alphabetical, each its own N-)
        private string BuildNSubstituentPrefix(List<string> subNames)
        {
            if (subNames.Count > 1 && subNames.Distinct().Count() == 1)
            {
                string locants = string.Join(",", Enumerable.Repeat("N", subNames.Count)); // "N,N"
                return locants + "-" + spec.numericalPrefixes[subNames.Count] + subNames[0]; // "N,N-dimethyl"
            }
            var sorted = subNames.OrderBy(w => w, StringComparer.Ordinal);
            return string.Join("-", sorted.Select(w => "N-" + w));                           // "N-ethyl-N-methyl"
        }

        private string[] NameSymmetricBridge(
            List<int> sideA, List<int> sideB,
            string substituentSuffix, string baseSuffix, string namePrefix)
        {
            SideInfo a = ClassifySide(sideA);
            SideInfo b = ClassifySide(sideB);
            if (a.Kind != SideKind.PlainAlkyl || b.Kind != SideKind.PlainAlkyl)
                throw new Exception("Ethers/amines with a carbonyl or ring side are not supported here");

            // Convention: shorter chain = substituent, longer = base chain
            int substituentLen = Math.Min(a.CarbonCount, b.CarbonCount);
            int baseLen = Math.Max(a.CarbonCount, b.CarbonCount);

            // namePrefix goes on AFTER FormatName, which lowercases the whole string —
            // this is what keeps the amine's "N-" uppercase.
            string body = FormatName(spec.alkylNames[substituentLen] + substituentSuffix
                                    + spec.alkylNames[baseLen] + baseSuffix);
            return new[] { namePrefix + body };
        }
        private string[] NameEther(List<int> sideA, List<int> sideB) =>
            NameSymmetricBridge(sideA, sideB, "oxy", spec.middle[""].name, "");
            // methoxymethane, methoxyethane, ethoxypropane — spec.middle[""] ("an|e")
            // rather than a hardcoded suffix

        private string[] NameSecondaryAmine(List<int> sideA, List<int> sideB) =>
            NameSymmetricBridge(sideA, sideB, "yl", "an|amine", "N-");
            // N-methylethanamine, N-ethylpropanamine

        private string[] NameEster(List<int> sideA, List<int> sideB)
        {
            SideInfo a = ClassifySide(sideA);
            SideInfo b = ClassifySide(sideB);
            var (acid, alkyl) = a.Kind == SideKind.AcidSide ? (a, b) : (b, a); // dispatch ⇒ exactly one AcidSide
            return new[] { AlkylSideName(alkyl) + " " + AcidSideName(acid, AcidClass.Ester) };
            // "methyl ethanoate", "ethyl propanoate" — methanoate (acidLen == 1) also works
        }

        private string[] NameAnhydride(List<int> sideA, List<int> sideB)
        {
            SideInfo a = ClassifySide(sideA); // dispatch ⇒ both AcidSide
            SideInfo b = ClassifySide(sideB);
            string stemA = AcidSideName(a, AcidClass.Anhydride); // "ethanoic"
            string stemB = AcidSideName(b, AcidClass.Anhydride);

            if (a.CarbonCount == b.CarbonCount)                  // symmetric ⇒ genuinely identical
                return new[] { stemA + " anhydride" };

            // Mixed: both acid names, alphabetical, then "anhydride".
            var ordered = new[] { stemA, stemB }.OrderBy(s => s, StringComparer.Ordinal).ToArray();
            return new[] { ordered[0] + " " + ordered[1] + " anhydride" };
        }

        private string[] NameSecondaryAmide(List<int> sideA, List<int> sideB)
        {
            SideInfo a = ClassifySide(sideA);
            SideInfo b = ClassifySide(sideB);
            var (acid, alkyl) = a.Kind == SideKind.AcidSide ? (a, b) : (b, a);

            string nPrefix = BuildNSubstituentPrefix(new List<string> { AlkylSideName(alkyl) }); // "N-methyl"
            return new[] { nPrefix + AcidSideName(acid, AcidClass.Amide) };                       // "N-methylethanamide"
        }

        private string[] NameTertiaryAmide(int bridgeIndex)
        {
            List<List<int>> sides = atoms.SplitAtBridgingAtomMultiway(bridgeIndex);
            List<SideInfo> infos = sides.Select(ClassifySide).ToList();

            SideInfo acid = infos.First(s => s.Kind == SideKind.AcidSide);          // dispatch ⇒ exactly one
            var nSubNames = infos.Where(s => s.Kind != SideKind.AcidSide)
                                 .Select(AlkylSideName).ToList();                    // the two N-substituents

            string nPrefix = BuildNSubstituentPrefix(nSubNames);                     // "N,N-dimethyl"
            return new[] { nPrefix + AcidSideName(acid, AcidClass.Amide) };          // "N,N-dimethylethanamide"
        }

        // ── Wave 2 / G1: unified substituent namer ─────────────────────────────────
        // One entry point for naming a substituent hanging off a parent skeleton (a
        // ring carbon or a chain junction). Returns a fully FormatName'd prefix word
        // ("methyl", "chloro", "(hydroxymethyl)", "(1-methylethyl)"), or throws — the
        // capability gate (Constraint 2: name it exactly or reject, never guess).
        // allowSuffixCapableTip: false on the aromatic path (D5) — a suffix-capable
        // tip there belongs to the parent/pattern tables, never to a prefix.
        private string NameSubstituent(int parentIndex, int attachIndex, HashSet<int> blocked,
                                       bool allowSuffixCapableTip)
        {
            if (atoms.BondOrder(parentIndex, attachIndex) != 1)
                throw new Exception("Substituents attached by a multiple bond are not supported");

            if (atoms.Atoms[attachIndex].Name != "Carbon")
                return NameHeteroatomSubstituent(attachIndex, blocked);

            return NameCarbonSubstituent(attachIndex, blocked, allowSuffixCapableTip);
        }

        // Degree-aware heteroatom naming — the anisole / N-methylphenylamine backstop
        // demanded by Further-Extension-Plan G1. In Wave 2 a non-terminal heteroatom
        // substituent is unreachable (bridging atoms and heteroatom–heteroatom pairs are
        // rejected upstream), so the degree guard is defence in depth; E6 replaces the
        // throw with alkoxy naming for O, and Wave 3 adds a PolyatomicGroup branch ABOVE
        // the degree guard for nitro.
        private string NameHeteroatomSubstituent(int attachIndex, HashSet<int> blocked)
        {
            if (atoms.AdjacentAtoms(attachIndex).Any(n => !blocked.Contains(n)))
                throw new Exception("Substituents extending beyond a single heteroatom are not supported");

            string formula = "C-" + atoms.Atoms[attachIndex].Symbol;
            if (spec.prefixOnly.ContainsKey(formula))
                return FormatName(spec.prefixOnly[formula]);
            if (spec.prefixOrSuffix.ContainsKey(formula))
                return FormatName(spec.prefixOrSuffix[formula].prefix);
            // Replaces GetSubstituentName's symbol.ToLower() fallback, which silently
            // emitted "cl" under specs with no halo entries (Q1). Reject instead.
            throw new Exception($"No prefix name available for substituent '{formula}'");
        }

        // Carbon substituent grammar (Wave 2, per ruling D4):
        //   spine = longest carbon path from the attachment atom (locant 1 = attachment),
        //   plus EITHER one tip functional group (on the last spine carbon only)
        //        OR bare methyl branches on spine carbons — never both.
        private string NameCarbonSubstituent(int attachIndex, HashSet<int> blocked,
                                             bool allowSuffixCapableTip)
        {
            List<int> subAtoms = atoms.CollectReachable(attachIndex, blocked); // subAtoms[0] == attachIndex

            // 1. Attachment guards: the substituent touches the parent skeleton exactly
            //    once, through attachIndex (rejects fused/bridged ring systems).
            if (atoms.AdjacentAtoms(attachIndex).Count(n => blocked.Contains(n)) != 1)
                throw new Exception("Fused or bridged ring systems are not supported");
            foreach (int a in subAtoms)
                if (a != attachIndex && atoms.AdjacentAtoms(a).Any(n => blocked.Contains(n)))
                    throw new Exception("Fused or bridged ring systems are not supported");

            // 2. Cycle guard (cyclohexyl / biphenyl substituents): the subgraph must be
            //    a tree. MUST run before LongestPathFrom, which assumes acyclicity.
            HashSet<int> subSet = new HashSet<int>(subAtoms);
            int directedEdges = subAtoms.Sum(a => atoms.AdjacentAtoms(a).Count(n => subSet.Contains(n)));
            if (directedEdges != 2 * (subAtoms.Count - 1))
                throw new Exception("Cyclic substituents are not supported");

            // 3. No unsaturation anywhere in the substituent. (The styrene shape is
            //    handled by the whole-molecule pattern table BEFORE G1 is called.)
            if (groups.Any(g => g is CarbonCarbonGroup && subAtoms.Any(g.Involves)))
                throw new Exception("Unsaturated substituents are not supported");

            // 4. Spine and inventory.
            HashSet<int> subCarbons = new HashSet<int>(
                subAtoms.Where(i => atoms.Atoms[i].Name == "Carbon"));
            List<int> spine = LongestPathFrom(attachIndex, subCarbons);
            int tip = spine[spine.Count - 1];

            // Heteroatoms may only sit on the spine tip.
            foreach (int a in subAtoms)
                if (atoms.Atoms[a].Name != "Carbon" && !atoms.AdjacentAtoms(tip).Contains(a))
                    throw new Exception("Substituents with functional groups on a non-tip atom are not supported");

            List<FunctionalGroup> tipGroups = groups
                .Where(g => !(g is CarbonCarbonGroup) && g.MainIndex == tip).ToList();
            List<int> branchCarbons = subCarbons.Where(c => !spine.Contains(c)).ToList();

            if (branchCarbons.Count == 0 && tipGroups.Count == 0)
                return FormatName(spec.alkylNames[spine.Count] + "yl");        // plain alkyl

            if (branchCarbons.Count == 0)                                      // tip-group alkyl
            {
                if (tipGroups.Count > 1)
                    throw new Exception("Substituents with multiple functional groups are not supported"); // Q2: -CF3, -CHCl2
                FunctionalGroup tg = tipGroups[0];
                // Carbon-subsuming groups (merged COOH/COCl/CON, and C≡N): their prefix
                // names include the tip carbon itself, so prefix+alkyl naming here is
                // structurally wrong (Q3). Reject; the chain path's candidate filter
                // then lets a suffix-bearing chain name the molecule correctly.
                if (tg is MergedGroup || tg.GroupFormula == "C≡N")
                    throw new Exception("Substituents containing carbon-based functional groups are not supported");
                string groupPrefix = TipGroupPrefix(tg.GroupFormula, allowSuffixCapableTip);
                return "(" + FormatName(groupPrefix + spec.alkylNames[spine.Count] + "yl") + ")";
            }

            // Branched alkyl: pure carbon, methyl branches only (D4).
            if (tipGroups.Count > 0 || subCarbons.Count != subAtoms.Count)
                throw new Exception("Branched substituents carrying functional groups are not supported");
            List<int> locants = new List<int>();
            foreach (int b in branchCarbons)
            {
                if (atoms.AlkylCounter(b) != 1)
                    throw new Exception("Substituent branches longer than methyl are not supported");
                int spineIndex = spine.FindIndex(s => atoms.AdjacentAtoms(b).Contains(s));
                if (spineIndex < 0)
                    throw new Exception("Substituent branches longer than methyl are not supported");
                locants.Add(spineIndex + 1);
            }
            string branchPart = NameSegment(new Dictionary<string, List<int>>
                { { spec.alkylNames[1] + "yl", locants } });                   // "1meth|ayl" / "1,1dimeth|ayl"
            return "(" + FormatName(branchPart + spec.alkylNames[spine.Count] + "yl") + ")";
        }

        // Prefix vocabulary for a (non-carbon-subsuming) tip functional group.
        private string TipGroupPrefix(string formula, bool allowSuffixCapable)
        {
            if (spec.prefixOnly.ContainsKey(formula))
                return spec.prefixOnly[formula];                               // chloro, bromo, ...
            if (spec.prefixOrSuffix.ContainsKey(formula))
            {
                if (!allowSuffixCapable)                                       // aromatic path (D5)
                    throw new Exception("Aromatic ring substituents carrying a principal-group tip are not supported");
                return spec.prefixOrSuffix[formula].prefix;                    // hydroxy, amino, ...
            }
            // e.g. tip C=O (endDependent only): was silent "" before (Q4). Reject.
            throw new Exception($"No prefix name available for group '{formula}'");
        }

        // Longest simple path from `start` through the given carbon set. The caller has
        // already verified the subgraph is a tree, so plain DFS terminates. Ties resolve
        // by adjacency order — deterministic, and equivalent under the methyl-only
        // branch rule (any leftover longer than methyl throws regardless of tie choice).
        private List<int> LongestPathFrom(int start, HashSet<int> allowed)
        {
            List<int> best = new List<int> { start };
            foreach (int n in atoms.AdjacentAtoms(start).Where(allowed.Contains))
            {
                HashSet<int> narrowed = new HashSet<int>(allowed);
                narrowed.Remove(start);
                List<int> tail = LongestPathFrom(n, narrowed);
                if (tail.Count + 1 > best.Count)
                {
                    best = new List<int> { start };
                    best.AddRange(tail);
                }
            }
            return best;
        }

        // ── Wave 2 / G4: aromatic substituent classification ──────────────────────
        private sealed class AromaticSubstituent
        {
            public int RingPosition;
            public string? ParentFormula;   // set for retained-parent candidates
            public string? ParentName;      // "phenol", "benzoic acid", ...
            public int Priority;            // spec suffix priority (parent candidates only)
            public string? ChainParentName; // "phenylmethanol" / "phenylethanone" / "ethenylbenzene"
            public string? PrefixName;      // G1 name (simple substituents only)
        }

        private AromaticSubstituent ClassifyAromaticSubstituent(
            List<int> ring, (int ringPosition, int atomIndex, bool isCarbon) sub, HashSet<int> ringSet)
        {
            int ringCarbon = ring[sub.ringPosition];

            if (!sub.isCarbon)
            {
                // Terminal, single-bonded heteroatom with a retained parent (O → phenol,
                // N → phenylamine). Candidacy also requires the live spec to rank the
                // group — under the Hydrocarbons spec nothing ranks, so the input falls
                // through to G1, which rejects cleanly.
                if (atoms.BondOrder(ringCarbon, sub.atomIndex) == 1
                    && !atoms.AdjacentAtoms(sub.atomIndex).Any(n => !ringSet.Contains(n)))
                {
                    string formula = "C-" + atoms.Atoms[sub.atomIndex].Symbol;
                    if (SpecificationData.AromaticHeteroatomParentNames.ContainsKey(formula)
                        && PriorityForFormula(formula) >= 0)
                        return new AromaticSubstituent
                        {
                            RingPosition = sub.ringPosition,
                            ParentFormula = formula,
                            ParentName = SpecificationData.AromaticHeteroatomParentNames[formula],
                            Priority = PriorityForFormula(formula)
                        };
                }
                return Simple();
            }

            List<int> subAtoms = atoms.CollectReachable(sub.atomIndex, ringSet);

            string? carbonParent = MatchCarbonParentFormula(sub.atomIndex, subAtoms);
            if (carbonParent != null && PriorityForFormula(carbonParent) >= 0)
                return new AromaticSubstituent
                {
                    RingPosition = sub.ringPosition,
                    ParentFormula = carbonParent,
                    ParentName = SpecificationData.AromaticCarbonParentNames[carbonParent],
                    Priority = PriorityForFormula(carbonParent)
                };

            string? chainParent = MatchChainParentPattern(sub.atomIndex, subAtoms);
            if (chainParent != null)
                return new AromaticSubstituent
                { RingPosition = sub.ringPosition, ChainParentName = chainParent };

            return Simple();

            AromaticSubstituent Simple() => new AromaticSubstituent
            {
                RingPosition = sub.ringPosition,
                PrefixName = NameSubstituent(ringCarbon, sub.atomIndex, ringSet,
                                             allowSuffixCapableTip: false) // throws if unnameable
            };
        }

        // A carbon parent is a ONE-carbon substituent whose only (non-CC) group is one
        // of the carbon-parent formulas, post-merge. Valence makes the shapes exact:
        // a COOH carbon has no spare bond, and a CHO carbon carrying Cl instead of H
        // has already merged into COCl.
        private string? MatchCarbonParentFormula(int attachIndex, List<int> subAtoms)
        {
            if (subAtoms.Count(i => atoms.Atoms[i].Name == "Carbon") != 1)
                return null;
            var attachGroups = groups
                .Where(g => !(g is CarbonCarbonGroup) && g.MainIndex == attachIndex).ToList();
            if (attachGroups.Count != 1)
                return null;
            string formula = attachGroups[0].GroupFormula;
            return SpecificationData.AromaticCarbonParentNames.ContainsKey(formula) ? formula : null;
        }

        // G4.4: whole-molecule chain-parent patterns. Exact-shape matches only —
        // anything similar-but-not-equal falls through to G1, which throws.
        private string? MatchChainParentPattern(int attachIndex, List<int> subAtoms)
        {
            int carbonCount = subAtoms.Count(i => atoms.Atoms[i].Name == "Carbon");
            var heteroGroups = groups
                .Where(g => !(g is CarbonCarbonGroup) && subAtoms.Contains(g.MainIndex)).ToList();
            bool hasCC = groups.Any(g => g is CarbonCarbonGroup && subAtoms.Any(g.Involves));

            // ring–CH2OH → phenylmethanol
            if (carbonCount == 1 && subAtoms.Count == 2 && !hasCC
                && heteroGroups.Count == 1 && heteroGroups[0].GroupFormula == "C-O")
                return SpecificationData.AromaticChainParentNames["CH2OH"];

            // ring–C(=O)CH3 → phenylethanone (the O's valence is full, so the second
            // carbon can only be bonded to the attach carbon — no extra check needed)
            if (carbonCount == 2 && subAtoms.Count == 3 && !hasCC
                && heteroGroups.Count == 1 && heteroGroups[0].GroupFormula == "C=O"
                && heteroGroups[0].MainIndex == attachIndex)
                return SpecificationData.AromaticChainParentNames["COCH3"];

            // ring–CH=CH2 → ethenylbenzene
            if (carbonCount == 2 && subAtoms.Count == 2 && heteroGroups.Count == 0
                && groups.Any(g => g is CarbonCarbonGroup cc && cc.GroupFormula == "C=C"
                                  && subAtoms.Contains(cc.MainIndex) && subAtoms.Contains(cc.OtherCarbonIndex)))
                return SpecificationData.AromaticChainParentNames["CHCH2"];

            return null;
        }

        // Rank / demote helpers — single source of truth is the live spec.
        private int PriorityForFormula(string formula)
        {
            if (spec.prefixOrSuffix.ContainsKey(formula)) return spec.prefixOrSuffix[formula].priority;
            if (spec.endDependentPrefixOrSuffix.ContainsKey(formula)) return spec.endDependentPrefixOrSuffix[formula].end.priority; // C=O → 6
            return -1;
        }

        private string PrefixForFormula(string formula)
        {
            if (spec.prefixOnly.ContainsKey(formula)) return spec.prefixOnly[formula];
            if (spec.prefixOrSuffix.ContainsKey(formula)) return spec.prefixOrSuffix[formula].prefix;      // hydroxy, amino, carboxy
            if (spec.endDependentPrefixOrSuffix.ContainsKey(formula)) return spec.endDependentPrefixOrSuffix[formula].end.prefix; // formyl
            throw new Exception($"No prefix name available for group '{formula}'");
        }

        // ── Phase 1: Cyclic naming ─────────────────────────────────────────────────
        private string[] NameCyclicMolecule()
        {
            // Wave 2 (D2): merging here is safe — dispatch has already run, so no
            // bridging heteroatom (and hence no ester/amide carbon the merge could
            // corrupt) can be present — and necessary, so ring substituents like -COOH
            // appear as one merged group for the parent table and G1's tip lookup.
            atoms.MergeFunctionalGroups(spec.merging);

            List<int> ring = atoms.FindRing();

            if (atoms.IsAromatic(ring))
                return NameAromaticMolecule(ring);

            HashSet<int> ringSet = new HashSet<int>(ring);
            if (groups.Any(g => g is CarbonCarbonGroup cc
                && ringSet.Contains(cc.MainIndex) && ringSet.Contains(cc.OtherCarbonIndex)))
                throw new Exception("Non-aromatic rings containing double or triple bonds are not supported");

            string ringBaseName = "cyclo" + spec.alkylNames[ring.Count];
            string middleName = spec.middle[""].name; // "an|e" for alkane (from spec)

            var substituents = atoms.GetRingSubstituents(ring);
            if (substituents.Count == 0)
                return new[] { FormatName(ringBaseName + middleName) };

            // G1 names every substituent — validation now lives inside the namer.
            var named = substituents
                .Select(s => (s.ringPosition,
                              name: NameSubstituent(ring[s.ringPosition], s.atomIndex, ringSet,
                                                    allowSuffixCapableTip: true)))
                .ToList();

            if (named.Count == 1)
                return new[] { FormatName(named[0].name + ringBaseName + middleName) };

            (int[] numbering, List<string> subNames) = FindBestRingNumbering(ring, named);
            return new[] { FormatName(BuildRingPrefixName(numbering, subNames) + ringBaseName + middleName) };
        }
        // ── Phase 2: Aromatic naming ───────────────────────────────────────────────
        private string[] NameAromaticMolecule(List<int> ring)
        {
            HashSet<int> ringSet = new HashSet<int>(ring);
            var substituents = atoms.GetRingSubstituents(ring);

            if (substituents.Count == 0)
                return new[] { "benzene" };

            // Enforce scope limit: ≤2 substituents   // ← DELETE THIS BLOCK IN PHASE D (E3)
            if (substituents.Count > 2)
                throw new Exception("Benzene with more than 2 substituents is not supported");

            var named = substituents
                .Select(s => (s.ringPosition,
                              name: NameSubstituent(ring[s.ringPosition], s.atomIndex, ringSet,
                                                    allowSuffixCapableTip: false)))
                .ToList();

            if (named.Count == 1)
                return new[] { FormatName(named[0].name + "benzene") };
                // e.g. "chlorobenzene", "methylbenzene"

            // 2 substituents: number to give lowest locants
            (int[] numbering, List<string> subNames) = FindBestRingNumbering(ring, named);
            string prefixName = BuildRingPrefixName(numbering, subNames);
            return new[] { FormatName(prefixName + "benzene") };
        }
        private (int[] bestNumbering, List<string> subNames) FindBestRingNumbering(
            List<int> ring,
            List<(int ringPosition, string name)> substituents,
            int anchorPos = -1) // G4: ring position forced to locant 1 (parent anchor)
        {
            int n = ring.Count;
            int[]? bestLocants = null;
            List<string>? bestNames = null;
            int bestStart = 0;
            int bestDirection = 0;

            foreach (int start in anchorPos >= 0 ? new[] { anchorPos } : Enumerable.Range(0, n).ToArray())
            {
                for (int dir = 0; dir < 2; dir++) // 0 = forward, 1 = reverse
                {
                    var pairs = new List<(int locant, string name)>();
                    foreach (var (ringPos, name) in substituents)
                    {
                        int locant = dir == 0
                            ? ((ringPos - start + n) % n) + 1
                            : ((start - ringPos + n) % n) + 1;
                        pairs.Add((locant, name));
                    }
                    pairs.Sort((a, b) => a.locant != b.locant
                        ? a.locant.CompareTo(b.locant)
                        : string.CompareOrdinal(AlphaKey(a.name), AlphaKey(b.name)));
                    int[] locants = pairs.Select(p => p.locant).ToArray();
                    List<string> nameList = pairs.Select(p => p.name).ToList();

                    bool isBetter = false;
                    if (bestLocants == null)
                    {
                        isBetter = true;
                    }
                    else
                    {
                        for (int i = 0; i < locants.Length; i++)
                        {
                            if (locants[i] < bestLocants[i]) { isBetter = true; break; }
                            if (locants[i] > bestLocants[i]) break;
                        }
                        // Tiebreaker: alphabetical order of names at first difference
                        if (!isBetter && locants.SequenceEqual(bestLocants))
                        {
                            for (int i = 0; i < nameList.Count; i++)
                            {
                                int cmp = string.CompareOrdinal(AlphaKey(nameList[i]), AlphaKey(bestNames![i]));
                                if (cmp < 0) { isBetter = true; break; }
                                if (cmp > 0) break;
                            }
                        }
                    }

                    if (isBetter)
                    {
                        bestLocants = locants;
                        bestNames = nameList;
                        bestStart = start;
                        bestDirection = dir;
                    }
                }
            }

            int[] numbering = new int[substituents.Count];
            List<string> subNames = new List<string>();
            for (int i = 0; i < substituents.Count; i++)
            {
                int ringPos = substituents[i].ringPosition;
                numbering[i] = bestDirection == 0
                    ? ((ringPos - bestStart + n) % n) + 1
                    : ((bestStart - ringPos + n) % n) + 1;
                subNames.Add(substituents[i].name);
            }
            return (numbering, subNames);
        }
        private string BuildRingPrefixName(int[] numbering, List<string> subNames)
        {
            Dictionary<string, List<int>> nameToLocants = new Dictionary<string, List<int>>();
            for (int i = 0; i < numbering.Length; i++)
            {
                if (nameToLocants.ContainsKey(subNames[i]))
                    nameToLocants[subNames[i]].Add(numbering[i]);
                else
                    nameToLocants[subNames[i]] = new List<int> { numbering[i] };
            }
            return NameSegment(nameToLocants);
        }
        public static void DisplaySpecDebug(string fileName) //no purpose besides displaying the rules extracted from the specification file
        {
            namingSpec spec = LoadSpecification(fileName);
            Console.WriteLine("alkylNames");
            ConsoleColor c = Console.ForegroundColor;
            Console.ForegroundColor = ConsoleColor.DarkRed;
            Console.WriteLine("----------");
            Console.ForegroundColor = c;
            foreach (string s in spec.alkylNames)
            {
                Console.Write("'" + s + "', ");
            }
            Console.WriteLine();
            Title("numericalPrefixes");
            foreach (string s in spec.numericalPrefixes)
            {
                Console.Write("'" + s + "', ");
            }
            Console.WriteLine();
            Title("merging");
            foreach (var v in spec.merging)
            {
                foreach (string s in v.Key)
                {
                    Console.Write("'" + s + "', ");
                }
                Console.Write(" : '" + v.Value + "'");
            }
            Console.WriteLine();
            Title("prefixOnly");
            foreach (var v in spec.prefixOnly)
            {
                Console.WriteLine("'" + v.Key + "' : '" + v.Value + "'");
            }
            Title("middle");
            foreach (var v in spec.middle)
            {
                Console.WriteLine("'" + v.Key + "' : '" + v.Value.name + "', '" + v.Value.priority + "'");
            }
            Title("prefixOrSuffix");
            foreach (var v in spec.prefixOrSuffix)
            {
                Console.WriteLine("'" + v.Key + "' : '" + v.Value.prefix + "', '" + v.Value.suffix + "', '" + v.Value.priority + "'");
            }
            Title("endDependentPrefixOrSuffix");
            foreach (var v in spec.endDependentPrefixOrSuffix)
            {
                Console.WriteLine("'" + v.Key + "' (MIDDLE) : '" + v.Value.middle.prefix + "', '" + v.Value.middle.suffix + "', '" + v.Value.middle.priority + "'");
                Console.WriteLine("'" + v.Key + "' (END) : '" + v.Value.end.prefix + "', '" + v.Value.end.suffix + "', '" + v.Value.end.priority + "'");
            }

            void Title(string s)
            {
                ConsoleColor c = Console.ForegroundColor;
                Console.ForegroundColor = ConsoleColor.DarkRed;
                Console.WriteLine("----------");
                Console.ForegroundColor = c;
                Console.WriteLine(s);
                Console.ForegroundColor = ConsoleColor.DarkRed;
                Console.WriteLine("----------");
                Console.ForegroundColor = c;
            }
        }
        public static void SaveSpecification(string fileName, namingSpec spec)
        {
            fileName += ".ExamSpec";
            using (BinaryWriter writeFile = new BinaryWriter(File.Open(fileName, FileMode.Create)))
            {
                writeFile.Write(((short)spec.alkylNames.Length)); //short is less than 32,000. that is enough
                foreach (string name in spec.alkylNames)
                {
                    writeFile.Write(name);
                }
                writeFile.Write((short)spec.numericalPrefixes.Length);
                foreach (string pref in spec.numericalPrefixes)
                {
                    writeFile.Write(pref);
                }
                writeFile.Write(((short)spec.merging.Count));
                foreach (KeyValuePair<HashSet<string>, string> v in spec.merging)
                {
                    writeFile.Write(((short)v.Key.Count));
                    foreach (string s in v.Key)
                    {
                        writeFile.Write(s);
                    }
                    writeFile.Write(v.Value);
                }
                writeFile.Write(((short)spec.prefixOnly.Count));
                foreach (KeyValuePair<string, string> v in spec.prefixOnly)
                {
                    writeFile.Write(v.Key);
                    writeFile.Write(v.Value);
                }
                List<(string key, string name)> sortedMiddle = new List<(string key, string name)>();
                sortedMiddle = spec.middle.Select(kv => (kv.Key, kv.Value.name)).ToList();
                sortedMiddle = sortedMiddle.OrderBy(x => spec.middle[x.key].priority).ToList(); //encode the 'priority' attribute as the order of the entires
                writeFile.Write(((short)sortedMiddle.Count));
                for (int i = sortedMiddle.Count - 1; i >= 0; i--)
                {
                    writeFile.Write(sortedMiddle[i].key);
                    writeFile.Write(sortedMiddle[i].name);
                }
                List<(string key, bool middleOnly, bool endOnly, string prefix, string suffix)> sortedSuffixes = spec.prefixOrSuffix.Select(kv => (kv.Key, false, false, kv.Value.prefix, kv.Value.suffix)).ToList();
                sortedSuffixes.AddRange(spec.endDependentPrefixOrSuffix.Select(kv => (kv.Key, true, false, kv.Value.middle.prefix, kv.Value.middle.suffix)).ToList()); //boolean second and third terms represent 'middle only' and 'end only' respectively
                sortedSuffixes.AddRange(spec.endDependentPrefixOrSuffix.Select(kv => (kv.Key, false, true, kv.Value.end.prefix, kv.Value.end.suffix)).ToList());
                sortedSuffixes = sortedSuffixes.OrderBy(x => getPriorityForSuffixSort(x)).ToList();
                writeFile.Write(((short)sortedSuffixes.Count));
                for (int i = sortedSuffixes.Count - 1; i >= 0; i--)
                {
                    writeFile.Write(sortedSuffixes[i].key);
                    writeFile.Write(sortedSuffixes[i].middleOnly);
                    writeFile.Write(sortedSuffixes[i].endOnly);
                    writeFile.Write(sortedSuffixes[i].prefix);
                    writeFile.Write(sortedSuffixes[i].suffix);
                }
                writeFile.Close();

                int getPriorityForSuffixSort((string key, bool middleOnly, bool endOnly, string prefix, string suffix) x)
                {
                    if (x.middleOnly)
                    {
                        return spec.endDependentPrefixOrSuffix[x.key].middle.priority;
                    }
                    else if (x.endOnly)
                    {
                        return spec.endDependentPrefixOrSuffix[x.key].end.priority;
                    }
                    else
                    {
                        return spec.prefixOrSuffix[x.key].priority;
                    }
                } //function because different entries stored in different dictionaries
            }
        }
        private static namingSpec LoadSpecification(string fileName)
        {
            fileName += ".ExamSpec";
            namingSpec spec = new namingSpec();
            using (BinaryReader readFile = new BinaryReader(File.Open(fileName, FileMode.Open)))
            {
                spec.alkylNames = new string[readFile.ReadInt16()];
                for (int i = 0; i < spec.alkylNames.Length; i++)
                {
                    spec.alkylNames[i] = readFile.ReadString();
                }
                spec.numericalPrefixes = new string[readFile.ReadInt16()];
                for (int i = 0; i < spec.numericalPrefixes.Length; i++)
                {
                    spec.numericalPrefixes[i] = readFile.ReadString();
                }
                spec.merging = new Dictionary<HashSet<string>, string>();
                int mergingCount = readFile.ReadInt16();
                for (int i = 0; i < mergingCount; i++)
                {
                    HashSet<string> toMerge = new HashSet<string>();
                    int toMergeCount = readFile.ReadInt16();
                    for (int j = 0; j < toMergeCount; j++)
                    {
                        toMerge.Add(readFile.ReadString());
                    }
                    spec.merging.Add(toMerge, readFile.ReadString());
                }
                spec.prefixOnly = new Dictionary<string, string>();
                int prefixOnlyCount = readFile.ReadInt16();
                for (int i = 0; i < prefixOnlyCount; i++)
                {
                    string key = readFile.ReadString();
                    string name = readFile.ReadString();
                    spec.prefixOnly.Add(key, name);
                }
                spec.middle = new Dictionary<string, (string name, int priority)>();
                int middleOnlyCount = readFile.ReadInt16();
                for (int i = middleOnlyCount - 1; i >= 0; i--)
                {
                    string key = readFile.ReadString();
                    string name = readFile.ReadString();
                    spec.middle.Add(key, (name, i + 1));
                }
                spec.prefixOrSuffix = new Dictionary<string, (string prefix, string suffix, int priority)>();
                spec.endDependentPrefixOrSuffix = new Dictionary<string, ((string prefix, string suffix, int priority) middle, (string prefix, string suffix, int priority) end)>();
                int suffixCount = readFile.ReadInt16();
                for (int i = suffixCount - 1; i >= 0; i--)
                {
                    string key = readFile.ReadString();
                    bool middleOnly = readFile.ReadBoolean();
                    bool endOnly = readFile.ReadBoolean();
                    string prefix = readFile.ReadString();
                    string suffix = readFile.ReadString();
                    if (middleOnly)
                    {
                        if (spec.endDependentPrefixOrSuffix.ContainsKey(key)) //end version of this entry already recorded, now ading the middle version
                        {
                            ((string prefix, string suffix, int priority) middle, (string prefix, string suffix, int priority) end) v = spec.endDependentPrefixOrSuffix[key];
                            v.middle.prefix = prefix;
                            v.middle.suffix = suffix;
                            v.middle.priority = i + 1;
                            spec.endDependentPrefixOrSuffix[key] = v; //cannot mutate complex values in a dictionary without overwriting the entry completely
                        }
                        else //this entry does not exist- make it new for this key
                        {
                            spec.endDependentPrefixOrSuffix.Add(key, ((prefix, suffix, i + 1), ("", "", -1)));
                        }
                    }
                    else if (endOnly)
                    {
                        if (spec.endDependentPrefixOrSuffix.ContainsKey(key)) //middle version of this entry already recorded, now ading the end version
                        {
                            ((string prefix, string suffix, int priority) middle, (string prefix, string suffix, int priority) end) v = spec.endDependentPrefixOrSuffix[key];
                            v.end.prefix = prefix;
                            v.end.suffix = suffix;
                            v.end.priority = i + 1;
                            spec.endDependentPrefixOrSuffix[key] = v;
                        }
                        else //this entry does not exist- make it new for this key
                        {
                            spec.endDependentPrefixOrSuffix.Add(key, (("", "", -1), (prefix, suffix, i + 1)));
                        }
                    }
                    else
                    {
                        if (spec.prefixOrSuffix.ContainsKey(key)) //already recorded somehow, overwrite it
                        {
                            (string prefix, string suffix, int priority) v = spec.prefixOrSuffix[key];
                            v.prefix = prefix;
                            v.suffix = suffix;
                            v.priority = i + 1;
                            spec.prefixOrSuffix[key] = v;
                        }
                        else
                        {
                            spec.prefixOrSuffix.Add(key, (prefix, suffix, i + 1));
                        }
                    }
                }
                readFile.Close();
            }
            return spec;
        }
        private void CheckGroups(List<FunctionalGroup> groups) //this check happpens after the merging step, so this checks the final groups before thay are named
        {
            foreach (FunctionalGroup g in groups)
            {
                if (!spec.prefixOnly.ContainsKey(g.GroupFormula) && !spec.prefixOrSuffix.ContainsKey(g.GroupFormula) && !spec.endDependentPrefixOrSuffix.ContainsKey(g.GroupFormula) && !spec.middle.ContainsKey(g.GroupFormula))
                {
                    throw new Exception($"{g.GroupFormula} is unrecognised"); //if none of the dictionaries have entries relating to this group
                }
            }
        }
        private void FindHighestPrioritySuffix()
        {
            // Only consider groups on the main chain — exclude groups on branch atoms
            HashSet<int> chainAtoms = new HashSet<int>(
                allChainsAndBranches.SelectMany(cab => cab.chain));

            int maxPriority = -1;
            foreach (FunctionalGroup fg in groups)
            {
                if (!chainAtoms.Contains(fg.MainIndex))
                    continue;
                int currentPriority;
                if (spec.prefixOrSuffix.ContainsKey(fg.GroupFormula))
                {
                    currentPriority = spec.prefixOrSuffix[fg.GroupFormula].priority;
                    if (currentPriority > maxPriority)
                    {
                        maxPriority = currentPriority;
                        suffixIsMiddle = false;
                        suffixIsEnd = false;
                        suffixFormula = fg.GroupFormula;
                        suffixRoot = spec.prefixOrSuffix[fg.GroupFormula].suffix; //these variables are set assuming this is the highest priority suffix. if not, these are overwritten next iter.
                    }
                }
                else if (spec.endDependentPrefixOrSuffix.ContainsKey(fg.GroupFormula)) //the group is enddependent
                {

                    if (atoms.IsAnEnd(fg.MainIndex)) //the group is on an end, so the end based info will be used
                    {
                        currentPriority = spec.endDependentPrefixOrSuffix[fg.GroupFormula].end.priority; //note that the end priority different than the middle one
                        if (currentPriority > maxPriority)
                        {
                            maxPriority = currentPriority;
                            suffixIsMiddle = false;
                            suffixIsEnd = true;
                            suffixFormula = fg.GroupFormula;
                            suffixRoot = spec.endDependentPrefixOrSuffix[fg.GroupFormula].end.suffix;
                        }
                    }
                    else //the group is in the middle, so the middle based info will be used
                    {
                        currentPriority = spec.endDependentPrefixOrSuffix[fg.GroupFormula].middle.priority;
                        if (currentPriority > maxPriority)
                        {
                            maxPriority = currentPriority;
                            suffixIsMiddle = true;
                            suffixIsEnd = false;
                            suffixFormula = fg.GroupFormula;
                            suffixRoot = spec.endDependentPrefixOrSuffix[fg.GroupFormula].middle.suffix;
                        }
                    }
                }
            }
        }
        private string FindHighestPriorityMiddleFormula()
        {
            int priority = -1;
            foreach (FunctionalGroup fg in groups)
            {
                if (spec.middle.ContainsKey(fg.GroupFormula))
                {
                    int currentPriority = spec.middle[fg.GroupFormula].priority;
                    if (currentPriority > priority)
                    {
                        priority = currentPriority;
                    }
                }
            }
            string formula = "";
            foreach (KeyValuePair<string, (string, int)> v in spec.middle)
            {
                if (v.Value.Item2 == priority)
                {
                    formula = v.Key;
                }
            }
            return formula;
        }
        public static void RestoreDefaultSpecifications() //only place with hard coded chemistry words
        {
            namingSpec fullSpec = new namingSpec();
            fullSpec.merging = new Dictionary<HashSet<string>, string>();
            fullSpec.prefixOnly = new Dictionary<string, string>();
            fullSpec.middle = new Dictionary<string, (string, int)>();
            fullSpec.prefixOrSuffix = new Dictionary<string, (string, string, int)>();
            fullSpec.endDependentPrefixOrSuffix = new Dictionary<string, ((string, string, int), (string, string, int))>();
            fullSpec.alkylNames = new string[] { "", "meth|a", "eth|a", "prop|a", "but|a", "pent|a", "hex|a", "hept|a", "oct|a", "non|a", "dec|a" }; //symbol to denote no double vowel is '|'
            fullSpec.numericalPrefixes = new string[] { "", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona" };
            fullSpec.merging.Add(new HashSet<string> { "C-O", "C=O" }, "COOH");
            fullSpec.merging.Add(new HashSet<string> { "C-Cl", "C=O" }, "COCl");
            fullSpec.merging.Add(new HashSet<string> { "C-N", "C=O" }, "CON");
            fullSpec.prefixOnly.Add("C-F", "fluoro");
            fullSpec.prefixOnly.Add("C-Cl", "chloro");
            fullSpec.prefixOnly.Add("C-Br", "bromo");
            fullSpec.prefixOnly.Add("C-I", "iodo");
            fullSpec.middle.Add("C=C", ("en|e", 2));
            fullSpec.middle.Add("C≡C", ("yn|e", 1));
            fullSpec.middle.Add("", ("an|e", 0));
            fullSpec.prefixOrSuffix.Add("COOH", ("carboxy", "oic acid", 10));
            fullSpec.prefixOrSuffix.Add("COCl", ("chlorocarbonyl", "oyl chloride", 9));
            fullSpec.prefixOrSuffix.Add("CON", ("carbamoyl", "amide", 8));
            fullSpec.prefixOrSuffix.Add("C≡N", ("cyano", "nitrile", 7));
            fullSpec.endDependentPrefixOrSuffix.Add("C=O", (("oxo", "one", 5), ("formyl", "al", 6)));
            fullSpec.prefixOrSuffix.Add("C-O", ("hydroxy", "ol", 4));
            fullSpec.prefixOrSuffix.Add("C-S", ("sulfanyl", "thiol", 3));
            fullSpec.prefixOrSuffix.Add("C-N", ("amino", "amine", 2));
            fullSpec.prefixOrSuffix.Add("C=N", ("imino", "imine", 1));
            //using this resource https://iupac.org/wp-content/uploads/2021/06/Organic-Brief-Guide-brochure_v1.1_June2021.pdf
            SaveSpecification("AllGroups", fullSpec);

            namingSpec hcarbonSpec = new namingSpec();
            hcarbonSpec.merging = new Dictionary<HashSet<string>, string>();
            hcarbonSpec.prefixOnly = new Dictionary<string, string>();
            hcarbonSpec.middle = new Dictionary<string, (string, int)>();
            hcarbonSpec.prefixOrSuffix = new Dictionary<string, (string, string, int)>();
            hcarbonSpec.endDependentPrefixOrSuffix = new Dictionary<string, ((string, string, int), (string, string, int))>();
            hcarbonSpec.alkylNames = new string[] { "", "meth|a", "eth|a", "prop|a", "but|a", "pent|a", "hex|a", "hept|a", "oct|a", "non|a", "dec|a" };
            hcarbonSpec.numericalPrefixes = new string[] { "", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona" };
            hcarbonSpec.middle.Add("C=C", ("en|e", 2));
            hcarbonSpec.middle.Add("C≡C", ("yn|e", 1));
            hcarbonSpec.middle.Add("", ("an|e", 0));
            SaveSpecification("Hydrocarbons", hcarbonSpec);

            namingSpec alkaneSpec = new namingSpec();
            alkaneSpec.merging = new Dictionary<HashSet<string>, string>();
            alkaneSpec.prefixOnly = new Dictionary<string, string>();
            alkaneSpec.middle = new Dictionary<string, (string, int)>();
            alkaneSpec.prefixOrSuffix = new Dictionary<string, (string, string, int)>();
            alkaneSpec.endDependentPrefixOrSuffix = new Dictionary<string, ((string, string, int), (string, string, int))>();
            alkaneSpec.alkylNames = new string[] { "", "meth|a", "eth|a", "prop|a", "but|a", "pent|a", "hex|a", "hept|a", "oct|a", "non|a", "dec|a" };
            alkaneSpec.numericalPrefixes = new string[] { "", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona" };
            alkaneSpec.middle.Add("", ("an", 0));
            SaveSpecification("Alkanes", alkaneSpec);
        }
    }
}