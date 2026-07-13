# Example Molecules

This directory contains example JSON structures for testing the Organic Namer API.

## Running the API Locally

### Prerequisites

1. **.NET 10 SDK** - Download from [dotnet.microsoft.com](https://dotnet.microsoft.com/download)
   ```bash
   # Verify installation
   dotnet --version
   ```

2. **Azure Functions Core Tools v4** - Required to run Azure Functions locally
   ```bash
   # macOS (Homebrew)
   brew tap azure/functions
   brew install azure-functions-core-tools@4

   # Windows (npm)
   npm install -g azure-functions-core-tools@4

   # Verify installation
   func --version
   ```

### Starting the API

1. Open a terminal and navigate to the Functions project:
   ```bash
   cd src/OrganicNamer.Functions
   ```

2. Build the project:
   ```bash
   dotnet build
   ```

3. Start the local Azure Functions host:
   ```bash
   func start
   ```

4. You should see output similar to:
   ```
   Azure Functions Core Tools

   Functions:
       Health: [GET] http://localhost:7071/api/health
       NameMolecule: [POST] http://localhost:7071/api/name
   ```

The API is now running at `http://localhost:7071`.

### Testing the API

**Check API health:**
```bash
curl http://localhost:7071/api/health
```

**Name a molecule using an example file:**
```bash
# From the repository root directory
curl -X POST http://localhost:7071/api/name \
  -H "Content-Type: application/json" \
  -d @examples/ethane.json
```

**Name a molecule with inline JSON:**
```bash
curl -X POST http://localhost:7071/api/name \
  -H "Content-Type: application/json" \
  -d '{
    "atoms": [
      {"element": "C", "bonds": [{"to": 1, "order": 1}]},
      {"element": "C", "bonds": [{"to": 0, "order": 1}]}
    ]
  }'
```

### Troubleshooting

| Issue | Solution |
|-------|----------|
| `func: command not found` | Install Azure Functions Core Tools (see prerequisites) |
| Port 7071 already in use | Use `func start --port 7072` to use a different port |
| Build errors | Run `dotnet restore` then `dotnet build` |
| "No functions found" | Ensure you're in the `src/OrganicNamer.Functions` directory |

## Examples

| File | Molecule | Expected Name | Description |
|------|----------|---------------|-------------|
| `ethane.json` | C₂H₆ | ethane | Simple alkane with 2 carbons |
| `ethanol.json` | C₂H₅OH | ethan-1-ol | Alcohol functional group |
| `propanone.json` | CH₃COCH₃ | propan-2-one | Ketone (acetone) |
| `2-methylpropane.json` | (CH₃)₃CH | 2-methylpropane | Branched alkane (isobutane) |
| `propene.json` | C₃H₆ | prop-1-ene | Alkene with double bond |
| `ethanoic-acid.json` | CH₃COOH | ethan-1-oic acid | Carboxylic acid (acetic acid) |
| `2-chloropropane.json` | CH₃CHClCH₃ | 2-chloropropane | Halogenated compound |
| `cyclohexane.json` | C₆H₁₂ | cyclohexane | Cycloalkane ring |
| `methylcyclohexane.json` | C₆H₁₁CH₃ | methylcyclohexane | Substituted cycloalkane |
| `methoxymethane.json` | CH₃OCH₃ | methoxymethane | Ether (heavy atoms only — H filled implicitly) |
| `methyl-ethanoate.json` | CH₃COOCH₃ | methyl ethanoate | Ester (heavy atoms only — H filled implicitly) |

## JSON Format

Each file follows this structure:
```json
{
  "atoms": [
    {
      "element": "C",
      "bonds": [
        {"to": 1, "order": 1}
      ]
    }
  ],
  "specificationSet": "AllGroups"
}
```

- **atoms**: Array of atoms with their element symbol and bonds
- **bonds**: Array of bonds where `to` is the atom index and `order` is 1 (single), 2 (double), or 3 (triple)
- **specificationSet**: "AllGroups", "Hydrocarbons", or "Alkanes"

### Specification Sets

| Set | Description | Supported Groups |
|-----|-------------|------------------|
| `AllGroups` | Full IUPAC nomenclature (default) | Alcohols, aldehydes, ketones, carboxylic acids, halogens, amines, etc. |
| `Hydrocarbons` | Carbon and hydrogen only | Alkanes, alkenes, alkynes |
| `Alkanes` | Single bonds only | Alkanes (saturated hydrocarbons) |

## API Responses

### Successful Response (200 OK)

```json
{
  "names": ["ethanol"]
}
```

The `names` array may contain multiple names for molecules with equivalent naming options.

### Error Response (422 Unprocessable Entity)

```json
{
  "error": "ChemistryError",
  "message": "The molecule seems to have no ends"
}
```

Common errors:
- **"Request body with atoms is required"** - Empty or missing atoms array
- **"The molecule seems to have no ends"** - Cyclic structure detected (not supported)
- **"Unknown element"** - Unrecognised element symbol

## Creating New Examples

### Step-by-Step Guide

1. **Draw your molecule on paper** - Number each atom starting from 0

2. **Start with carbon atoms** - List carbons first (indices 0, 1, 2, ...)

3. **Add other atoms** - Hydrogens, oxygen, nitrogen, halogens, etc.

4. **Define bonds bidirectionally** - If atom 0 bonds to atom 1, both atoms must list the bond:
   ```json
   {"element": "C", "bonds": [{"to": 1, "order": 1}]},  // atom 0
   {"element": "C", "bonds": [{"to": 0, "order": 1}]}   // atom 1
   ```

5. **Use correct bond orders**:
   - `1` = single bond (C-C, C-H, C-O)
   - `2` = double bond (C=C, C=O)
   - `3` = triple bond (C≡C)

6. **Check valencies**:
   - Carbon: 4 bonds total
   - Hydrogen: 1 bond
   - Oxygen: 2 bonds
   - Nitrogen: 3 bonds
   - Halogens (Cl, Br, F, I): 1 bond

### Example: Building Methanol (CH₃OH)

```
Molecule:     H
              |
         H -- C -- O -- H
              |
              H

Atom indices: H(1)
              |
        H(2)--C(0)--O(3)--H(4)
              |
              H(5)
```

```json
{
  "atoms": [
    {"element": "C", "bonds": [{"to": 1, "order": 1}, {"to": 2, "order": 1}, {"to": 3, "order": 1}, {"to": 5, "order": 1}]},
    {"element": "H", "bonds": [{"to": 0, "order": 1}]},
    {"element": "H", "bonds": [{"to": 0, "order": 1}]},
    {"element": "O", "bonds": [{"to": 0, "order": 1}, {"to": 4, "order": 1}]},
    {"element": "H", "bonds": [{"to": 3, "order": 1}]},
    {"element": "H", "bonds": [{"to": 0, "order": 1}]}
  ],
  "specificationSet": "AllGroups"
}
```

### Validation Tips

- Count total bonds for each atom to verify valency
- Ensure every bond appears in both atoms' bond lists
- Test with the API before committing
