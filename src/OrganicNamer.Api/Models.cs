using OrganicNamer.Core;

namespace OrganicNamer.Api
{
    // Input model - uses ElementGraph.AtomInput from Core
    public class MoleculeRequest
    {
        public List<ElementGraph.AtomInput> Atoms { get; set; } = new();
        public string? SpecificationSet { get; set; } = "AllGroups";
    }

    // Output model
    public class MoleculeResponse
    {
        public List<string> Names { get; set; } = new();
    }

    // Error model
    public class ErrorResponse
    {
        public string Error { get; set; } = "";
        public string Message { get; set; } = "";
    }
}
