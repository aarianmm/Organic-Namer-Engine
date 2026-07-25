using System.Reflection;
using Microsoft.OpenApi;
using OrganicNamer.Core;
using OrganicNamer.Api;

var builder = WebApplication.CreateBuilder(args);

builder.Services.AddCors(options =>
{
    options.AddDefaultPolicy(policy =>
        policy.WithOrigins("https://thealkane.com", "https://www.thealkane.com")
            .AllowAnyMethod()
            .AllowAnyHeader());
});

builder.Services.AddEndpointsApiExplorer();
builder.Services.AddSwaggerGen(options =>
{
    options.SwaggerDoc("v1", new OpenApiInfo
    {
        Title = "OrganicNamer API",
        Version = "v1",
        Description = "Generates IUPAC names for organic molecules described as a graph of atoms and bonds."
    });

    foreach (var xmlFile in new[] { "OrganicNamer.Api.xml", "OrganicNamer.Core.xml" })
    {
        var xmlPath = Path.Combine(AppContext.BaseDirectory, xmlFile);
        if (File.Exists(xmlPath))
            options.IncludeXmlComments(xmlPath);
    }
});

var app = builder.Build();

app.UseCors();

app.UseSwagger();
app.UseSwaggerUI(options =>
{
    options.SwaggerEndpoint("/swagger/v1/swagger.json", "OrganicNamer API v1");
    options.RoutePrefix = "swagger";
});

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
})
.WithName("NameMolecule")
.WithTags("Naming")
.WithSummary("Generate the IUPAC name for a molecule")
.WithDescription("Accepts a molecule as a list of atoms with bonds to other atoms by index. " +
    "Hydrogens may be omitted and are filled in automatically to satisfy valency. " +
    "Returns one name, or several when the molecule has equally valid alternative names.")
.Produces<MoleculeResponse>(StatusCodes.Status200OK)
.Produces<ErrorResponse>(StatusCodes.Status400BadRequest)
.Produces<ErrorResponse>(StatusCodes.Status422UnprocessableEntity);

app.MapGet("/api/health", () => Results.Ok(new
{
    Status = "healthy",
    Version = "1.0.0",
    Specifications = new[] { "AllGroups", "Hydrocarbons", "Alkanes" }
}))
.WithName("Health")
.WithTags("Health")
.WithSummary("Report service health and supported specification sets");

app.Run();
