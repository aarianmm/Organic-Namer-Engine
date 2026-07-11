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
    }
}
