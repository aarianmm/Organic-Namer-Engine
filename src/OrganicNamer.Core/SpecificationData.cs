namespace OrganicNamer.Core
{
    public static class SpecificationData
    {
        // Periodic table - from RestoreDefaultPeriodicTable
        public static Dictionary<string, (string name, int valency)> PeriodicTable { get; } =
            new Dictionary<string, (string, int)>
            {
                { "C", ("Carbon", 4) },
                { "H", ("Hydrogen", 1) },
                { "O", ("Oxygen", 2) },
                { "N", ("Nitrogen", 3) },
                { "F", ("Fluorine", 1) },
                { "Cl", ("Chlorine", 1) },
                { "Br", ("Bromine", 1) },
                { "I", ("Iodine", 1) },
                { "S", ("Sulphur", 2) }
            };

        // Specifications - from RestoreDefaultSpecifications
        public static IUPAC.namingSpec AllGroups { get; } = BuildAllGroupsSpec();
        public static IUPAC.namingSpec Hydrocarbons { get; } = BuildHydrocarbonsSpec();
        public static IUPAC.namingSpec Alkanes { get; } = BuildAlkanesSpec();

        private static IUPAC.namingSpec BuildAllGroupsSpec()
        {
            var spec = new IUPAC.namingSpec();
            spec.merging = new Dictionary<HashSet<string>, string>();
            spec.prefixOnly = new Dictionary<string, string>();
            spec.middle = new Dictionary<string, (string, int)>();
            spec.prefixOrSuffix = new Dictionary<string, (string, string, int)>();
            spec.endDependentPrefixOrSuffix = new Dictionary<string, ((string, string, int), (string, string, int))>();

            spec.alkylNames = new string[] { "", "meth|a", "eth|a", "prop|a", "but|a", "pent|a", "hex|a", "hept|a", "oct|a", "non|a", "dec|a" };
            spec.numericalPrefixes = new string[] { "", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona" };

            spec.merging.Add(new HashSet<string> { "C-O", "C=O" }, "COOH");
            spec.merging.Add(new HashSet<string> { "C-Cl", "C=O" }, "COCl");
            spec.merging.Add(new HashSet<string> { "C-N", "C=O" }, "CON");

            spec.prefixOnly.Add("C-F", "fluoro");
            spec.prefixOnly.Add("C-Cl", "chloro");
            spec.prefixOnly.Add("C-Br", "bromo");
            spec.prefixOnly.Add("C-I", "iodo");
            spec.prefixOnly.Add("NO2", "nitro");   // Wave 3 / E4 (AllGroups only - the spec gate)

            spec.middle.Add("C=C", ("en|e", 2));
            spec.middle.Add("C≡C", ("yn|e", 1));
            spec.middle.Add("", ("an|e", 0));

            spec.prefixOrSuffix.Add("COOH", ("carboxy", "oic acid", 10));
            spec.prefixOrSuffix.Add("COCl", ("chlorocarbonyl", "oyl chloride", 9));
            spec.prefixOrSuffix.Add("CON", ("carbamoyl", "amide", 8));
            spec.prefixOrSuffix.Add("C≡N", ("cyano", "nitrile", 7));
            spec.endDependentPrefixOrSuffix.Add("C=O", (("oxo", "one", 5), ("formyl", "al", 6)));
            spec.prefixOrSuffix.Add("C-O", ("hydroxy", "ol", 4));
            spec.prefixOrSuffix.Add("C-S", ("sulfanyl", "thiol", 3));
            spec.prefixOrSuffix.Add("C-N", ("amino", "amine", 2));
            spec.prefixOrSuffix.Add("C=N", ("imino", "imine", 1));

            return spec;
        }

        private static IUPAC.namingSpec BuildHydrocarbonsSpec()
        {
            var spec = new IUPAC.namingSpec();
            spec.merging = new Dictionary<HashSet<string>, string>();
            spec.prefixOnly = new Dictionary<string, string>();
            spec.middle = new Dictionary<string, (string, int)>();
            spec.prefixOrSuffix = new Dictionary<string, (string, string, int)>();
            spec.endDependentPrefixOrSuffix = new Dictionary<string, ((string, string, int), (string, string, int))>();

            spec.alkylNames = new string[] { "", "meth|a", "eth|a", "prop|a", "but|a", "pent|a", "hex|a", "hept|a", "oct|a", "non|a", "dec|a" };
            spec.numericalPrefixes = new string[] { "", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona" };

            spec.middle.Add("C=C", ("en|e", 2));
            spec.middle.Add("C≡C", ("yn|e", 1));
            spec.middle.Add("", ("an|e", 0));

            return spec;
        }

        private static IUPAC.namingSpec BuildAlkanesSpec()
        {
            var spec = new IUPAC.namingSpec();
            spec.merging = new Dictionary<HashSet<string>, string>();
            spec.prefixOnly = new Dictionary<string, string>();
            spec.middle = new Dictionary<string, (string, int)>();
            spec.prefixOrSuffix = new Dictionary<string, (string, string, int)>();
            spec.endDependentPrefixOrSuffix = new Dictionary<string, ((string, string, int), (string, string, int))>();

            spec.alkylNames = new string[] { "", "meth|a", "eth|a", "prop|a", "but|a", "pent|a", "hex|a", "hept|a", "oct|a", "non|a", "dec|a" };
            spec.numericalPrefixes = new string[] { "", "", "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona" };

            spec.middle.Add("", ("an", 0));

            return spec;
        }

        // ── Wave 2 / G4: retained aromatic parent names (E1) ──────────────────────
        // Board-exact strings (phenylamine, NEVER aniline). Keyed by the substituent's
        // post-merge group formula. Split hetero-attached vs carbon-attached (D6) so a
        // C-O TIP group on a carbon substituent (-CH2OH) can never match phenol.
        // Priorities are NOT stored here — they come from the live namingSpec.
        public static IReadOnlyDictionary<string, string> AromaticHeteroatomParentNames { get; } =
            new Dictionary<string, string>
            {
                { "C-O", "phenol" },
                { "C-N", "phenylamine" },
            };

        public static IReadOnlyDictionary<string, string> AromaticCarbonParentNames { get; } =
            new Dictionary<string, string>
            {
                { "COOH", "benzoic acid" },
                { "COCl", "benzoyl chloride" },
                { "CON",  "benzamide" },       // primary aromatic amide — ring-substituent
                                               // pattern, owned here (secondary = E6)
                { "C≡N",  "benzonitrile" },
                { "C=O",  "benzaldehyde" },
            };

        // G4.4: whole-molecule chain-parent names — mono-substituted benzene only.
        public static IReadOnlyDictionary<string, string> AromaticChainParentNames { get; } =
            new Dictionary<string, string>
            {
                { "CH2OH", "phenylmethanol" },
                { "COCH3", "phenylethanone" },
                { "CHCH2", "ethenylbenzene" },
            };

        // ── Wave 4 / E6: retained strings for bridged-path ring sides ──────────────
        // ("phenylamine" is reused from AromaticHeteroatomParentNames["C-N"] above —
        // single source of truth, per G4.)
        public const string AromaticAlkylSideName = "phenyl";
        public const string AromaticEsterAcidStem = "benzoate";
        public const string AromaticAmideAcidStem = "benzamide";
    }
}
