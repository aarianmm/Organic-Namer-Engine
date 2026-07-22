FROM mcr.microsoft.com/dotnet/sdk:10.0 AS build
WORKDIR /src

COPY src/OrganicNamer.Core/OrganicNamer.Core.csproj src/OrganicNamer.Core/
COPY src/OrganicNamer.Api/OrganicNamer.Api.csproj src/OrganicNamer.Api/
RUN dotnet restore src/OrganicNamer.Api/OrganicNamer.Api.csproj

COPY src/ src/
RUN dotnet publish src/OrganicNamer.Api/OrganicNamer.Api.csproj -c Release -o /app --no-restore

FROM mcr.microsoft.com/dotnet/aspnet:10.0 AS runtime
WORKDIR /app
COPY --from=build /app .

ENV ASPNETCORE_HTTP_PORTS=8080
EXPOSE 8080
ENTRYPOINT ["dotnet", "OrganicNamer.Api.dll"]
