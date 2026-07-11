namespace OrganicNamer.Core
{
    public class Element
    {
        private int[] bondIndexes;
        public int[] BondIndexes { get { return bondIndexes; } }
        private string name;
        public string Name { get { return name; } }
        private string symbol;
        public string Symbol { get { return symbol; } }
        public int MaxBonds { get { return bondIndexes.Length; } }
        public int BondsFree { get { return Array.FindAll(bondIndexes, i => i == -1).Length; } }
        public Element(string symbol, string name, int valency)
        {
            this.symbol = symbol;
            this.name = name;
            bondIndexes = new int[valency];
            Array.Fill(bondIndexes, -1);
        }
        public void AddHalfBond(int bondIndex, int bondOrder) //bond order; 1 = single, 2 = double etc
        {
            for (int i = 0; i < bondOrder; i++)
            {
                bondIndexes[Array.FindIndex(bondIndexes, i => i == -1)] = bondIndex; //https://stackoverflow.com/questions/5400895/using-c-sharp-linq-to-return-first-index-of-null-empty-occurrence-in-an-array
            }
        }
    }
}