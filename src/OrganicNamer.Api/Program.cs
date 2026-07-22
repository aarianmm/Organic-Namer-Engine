using OrganicNamer.Core;
using OrganicNamer.Api;

var builder = WebApplication.CreateBuilder(args);

builder.Services.AddCors(options =>
{
    options.AddDefaultPolicy(policy =>
        policy.WithOrigins("https://alkane.uk", "https://www.alkane.uk")
            .AllowAnyMethod()
            .AllowAnyHeader());
});

var app = builder.Build();

app.UseCors();

app.MapPost("/api/name", (MoleculeRequest request) =>
{
    if (request.Atoms == null || request.Atoms.Count == 0)
        return Results.BadRequest(new ErrorResponse
        {
            Error = "ValidationError",
            Message = "Request body with atoms is required"
        });

    try
    {
        var spec = request.SpecificationSet switch
        {
            "Hydrocarbons" => SpecificationData.Hydrocarbons,
            "Alkanes"      => SpecificationData.Alkanes,
            _              => SpecificationData.AllGroups
        };

        var filledAtoms = ElementGraph.FillImplicitHydrogens(request.Atoms, SpecificationData.PeriodicTable);
        var graph = ElementGraph.FromJsonAtoms(filledAtoms, SpecificationData.PeriodicTable);
        var namer = new IUPAC(spec, graph);

        return Results.Ok(new MoleculeResponse { Names = namer.names.ToList() });
    }
    catch (Exception ex)
    {
        return Results.UnprocessableEntity(new ErrorResponse
        {
            Error = "ChemistryError",
            Message = ex.Message
        });
    }
});

app.MapGet("/api/health", () => Results.Ok(new
{
    Status = "healthy",
    Version = "1.0.0",
    Specifications = new[] { "AllGroups", "Hydrocarbons", "Alkanes" }
}));

app.Run();
