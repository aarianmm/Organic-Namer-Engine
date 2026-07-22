using OrganicNamer.Core;

namespace OrganicNamer.Api
{
    /// <summary>Request body for <c>POST /api/name</c>: a molecule described as a graph of atoms and bonds.</summary>
    public class MoleculeRequest
    {
        /// <summary>The atoms making up the molecule. Implicit hydrogens are filled in automatically, so they may be omitted.</summary>
        public List<ElementGraph.AtomInput> Atoms { get; set; } = new();

        /// <summary>Which functional groups the namer is allowed to recognize. One of "AllGroups", "Hydrocarbons", or "Alkanes". Defaults to "AllGroups".</summary>
        public string? SpecificationSet { get; set; } = "AllGroups";
    }

    /// <summary>Successful response from <c>POST /api/name</c>.</summary>
    public class MoleculeResponse
    {
        /// <summary>The generated IUPAC name(s). More than one entry means the molecule has equivalent valid names.</summary>
        public List<string> Names { get; set; } = new();
    }

    /// <summary>Error response returned with a non-2xx status code.</summary>
    public class ErrorResponse
    {
        /// <summary>Machine-readable error category, e.g. "ValidationError" or "ChemistryError".</summary>
        public string Error { get; set; } = "";

        /// <summary>Human-readable description of what went wrong.</summary>
        public string Message { get; set; } = "";
    }
}
