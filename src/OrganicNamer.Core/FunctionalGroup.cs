namespace OrganicNamer.Core
{
    public class FunctionalGroup
    {
        protected string elementSymbol;
        public string ElementSymbol { get { return elementSymbol; } }
        protected string bondSymbol;
        public virtual string GroupFormula { get { return "C" + bondSymbol + elementSymbol; } } //ie C=O or C≣N
        protected int mainIndex;
        public int MainIndex { get { return mainIndex; } }
        public FunctionalGroup(string elementSymbol, int bondOrder, int mainIndex)
        {
            this.elementSymbol = elementSymbol;
            this.mainIndex = mainIndex;
            switch (bondOrder)
            {
                case 1:
                    bondSymbol = "-";
                    break;
                case 2:
                    bondSymbol = "=";
                    break;
                case 3:
                    bondSymbol = "≡";
                    break;
                case 4:
                    bondSymbol = "≣"; //above 4 never used in real organic chemistry
                    break;
                default:
                    if (bondOrder > 0) //no symbol for 5+ lines, so the number is simply used
                    {
                        bondSymbol = bondOrder.ToString();
                    }
                    else
                    {
                        throw new Exception("Invalid bond order");
                    }                    
                    break;
            }
        }
        virtual public bool Involves(int index)
        {
            return mainIndex == index;
        }
    }
    public class CarbonCarbonGroup : FunctionalGroup
    {
        int otherCarbonIndex;
        public int OtherCarbonIndex { get { return otherCarbonIndex; } }
        public CarbonCarbonGroup(int bondOrder, int carbonIndex, int otherCarbonIndex) : base("C", bondOrder, carbonIndex)
        {
            this.otherCarbonIndex = otherCarbonIndex;
        }
        public bool IsSame(CarbonCarbonGroup otherGroup)
        {
            return (otherGroup.mainIndex == mainIndex && otherGroup.otherCarbonIndex == otherCarbonIndex) || (otherGroup.mainIndex == otherCarbonIndex && otherGroup.otherCarbonIndex == mainIndex);
        }
        override public bool Involves(int index)
        {
            return mainIndex == index || otherCarbonIndex == index;
        }
    }
    public class MergedGroup : FunctionalGroup
    {
        string formula;
        public override string GroupFormula { get { return formula; } }
        public MergedGroup(string formula, int carbonIndex) : base("", 1, carbonIndex)
        {
            this.formula = formula;
        }
    }
    public class CarbonOtherCarbonGroup : FunctionalGroup //future proofing
    {
        int[] carbonIndexes;
        public override string GroupFormula { get { return "C" + bondSymbol + elementSymbol + bondSymbol + "C"; } }
        public CarbonOtherCarbonGroup(string centralElementSymbol, int centralAtomIndex, int[] carbonIndexes) : base(centralElementSymbol, 1, centralAtomIndex)
        {
            this.carbonIndexes = carbonIndexes;
        }
        override public bool Involves(int index)
        {
            return carbonIndexes.Contains(index);
        }
    }
    public class PolyatomicGroup : FunctionalGroup //G2 (Wave 3): multi-heteroatom pattern, e.g. NO2
    {
        int[] memberIndexes; //the pattern's heteroatoms (nitro: N, O, O) - NOT the attached carbon
        string formula;
        public override string GroupFormula { get { return formula; } }
        public IReadOnlyList<int> MemberIndexes { get { return memberIndexes; } }
        public PolyatomicGroup(string formula, int carbonIndex, int[] memberIndexes) : base("", 1, carbonIndex)
        {
            this.formula = formula;
            this.memberIndexes = memberIndexes;
        }
        override public bool Involves(int index)
        {
            return mainIndex == index || memberIndexes.Contains(index);
        }
    }
}