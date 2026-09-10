using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Text.Json;
using System.Text.Json.Nodes;
using OrganicNamer.Core;

namespace OrganicNamer.Benchmarks;

internal static class Program
{
    private static readonly JsonSerializerOptions AtomJsonOptions = new()
    {
        PropertyNameCaseInsensitive = true
    };

    private static int Main(string[] args)
    {
        if (args.Length == 0)
        {
            PrintUsage();
            return 1;
        }

        var mode = args[0];
        var options = ParseOptions(args.Skip(1).ToArray());

        try
        {
            switch (mode)
            {
                case "name":
                    RunNameMode(options);
                    return 0;
                case "bench":
                    RunBenchMode(options);
                    return 0;
                default:
                    Console.Error.WriteLine($"Unknown mode '{mode}'. Expected 'name' or 'bench'.");
                    PrintUsage();
                    return 1;
            }
        }
        catch (RunnerUsageException ex)
        {
            Console.Error.WriteLine(ex.Message);
            PrintUsage();
            return 1;
        }
    }

    private static void PrintUsage()
    {
        Console.Error.WriteLine("Usage:");
        Console.Error.WriteLine("  dotnet run -c Release --project benchmarks/runner/OrganicNamer.Benchmarks -- name --in <input.jsonl> --out <output.jsonl>");
        Console.Error.WriteLine("  dotnet run -c Release --project benchmarks/runner/OrganicNamer.Benchmarks -- bench --in <input.jsonl> --out <results.json> [--warmup 3] [--reps 10]");
    }

    private sealed class RunnerUsageException(string message) : Exception(message);

    private sealed class Options
    {
        public string? In;
        public string? Out;
        public int Warmup = 3;
        public int Reps = 10;
    }

    private static Options ParseOptions(string[] args)
    {
        var options = new Options();
        for (int i = 0; i < args.Length; i++)
        {
            switch (args[i])
            {
                case "--in":
                    options.In = RequireValue(args, ref i, "--in");
                    break;
                case "--out":
                    options.Out = RequireValue(args, ref i, "--out");
                    break;
                case "--warmup":
                    options.Warmup = int.Parse(RequireValue(args, ref i, "--warmup"));
                    break;
                case "--reps":
                    options.Reps = int.Parse(RequireValue(args, ref i, "--reps"));
                    break;
                default:
                    throw new RunnerUsageException($"Unknown argument '{args[i]}'.");
            }
        }

        if (options.In == null) throw new RunnerUsageException("--in is required.");
        if (options.Out == null) throw new RunnerUsageException("--out is required.");
        return options;
    }

    private static string RequireValue(string[] args, ref int i, string flag)
    {
        if (i + 1 >= args.Length) throw new RunnerUsageException($"{flag} requires a value.");
        i++;
        return args[i];
    }

    // ------------------------------------------------------------------
    // Shared molecule parsing
    // ------------------------------------------------------------------

    private enum MoleculeStatus { Named, Rejected, Invalid }

    private sealed class ParsedMolecule
    {
        public required JsonObject Passthrough;
        public List<ElementGraph.AtomInput>? Atoms;
        public IUPAC.namingSpec? Spec;
        public MoleculeStatus? PreValidationStatus; // set if invalid before we even try naming
        public string? PreValidationError;
    }

    // Parses a single JSONL line into a ParsedMolecule. Throws only for genuinely
    // unparsable JSON (caller turns that into an "invalid" record too).
    private static ParsedMolecule ParseLine(string line)
    {
        JsonNode? node;
        try
        {
            node = JsonNode.Parse(line);
        }
        catch (Exception ex)
        {
            return new ParsedMolecule
            {
                Passthrough = new JsonObject(),
                PreValidationStatus = MoleculeStatus.Invalid,
                PreValidationError = $"Malformed JSON: {ex.Message}"
            };
        }

        if (node is not JsonObject obj)
        {
            return new ParsedMolecule
            {
                Passthrough = new JsonObject(),
                PreValidationStatus = MoleculeStatus.Invalid,
                PreValidationError = "Line is not a JSON object."
            };
        }

        // Clone so we can freely hand the original node tree back out untouched.
        var passthrough = obj.DeepClone().AsObject();

        List<ElementGraph.AtomInput>? atoms = null;
        string? atomsError = null;
        if (!obj.TryGetPropertyValue("atoms", out var atomsNode) || atomsNode == null)
        {
            atomsError = "Missing required field 'atoms'.";
        }
        else
        {
            try
            {
                atoms = atomsNode.Deserialize<List<ElementGraph.AtomInput>>(AtomJsonOptions);
            }
            catch (Exception ex)
            {
                atomsError = $"Could not parse 'atoms': {ex.Message}";
            }
        }

        if (atomsError == null && (atoms == null || atoms.Count == 0))
        {
            atomsError = "'atoms' must be a non-empty array.";
        }

        if (atomsError == null && atoms != null)
        {
            for (int i = 0; i < atoms.Count; i++)
            {
                foreach (var bond in atoms[i].Bonds)
                {
                    if (bond.To < 0 || bond.To >= atoms.Count)
                    {
                        atomsError = $"Atom {i} has a bond to out-of-range index {bond.To} (atom count = {atoms.Count}).";
                        break;
                    }
                }
                if (atomsError != null) break;
            }
        }

        if (atomsError != null)
        {
            return new ParsedMolecule
            {
                Passthrough = passthrough,
                PreValidationStatus = MoleculeStatus.Invalid,
                PreValidationError = atomsError
            };
        }

        string specName = "AllGroups";
        if (obj.TryGetPropertyValue("specificationSet", out var specNode) && specNode != null)
        {
            specName = specNode.GetValue<string>();
        }

        var spec = specName switch
        {
            "Hydrocarbons" => SpecificationData.Hydrocarbons,
            "Alkanes" => SpecificationData.Alkanes,
            _ => SpecificationData.AllGroups
        };

        return new ParsedMolecule
        {
            Passthrough = passthrough,
            Atoms = atoms,
            Spec = spec
        };
    }

    // Runs the actual naming pipeline (the same 3-line sequence OrganicNamer.Api uses),
    // timing only that pipeline. Never throws — chemistry-engine exceptions are caught
    // and reported as a rejection.
    private static (MoleculeStatus status, string[] names, string? error, double micros) NameMolecule(
        List<ElementGraph.AtomInput> atoms, IUPAC.namingSpec spec)
    {
        long start = Stopwatch.GetTimestamp();
        try
        {
            var filledAtoms = ElementGraph.FillImplicitHydrogens(atoms, SpecificationData.PeriodicTable);
            var graph = ElementGraph.FromJsonAtoms(filledAtoms, SpecificationData.PeriodicTable);
            var namer = new IUPAC(spec, graph);
            long end = Stopwatch.GetTimestamp();
            double micros = Stopwatch.GetElapsedTime(start, end).TotalMicroseconds;
            return (MoleculeStatus.Named, namer.names, null, micros);
        }
        catch (Exception ex)
        {
            long end = Stopwatch.GetTimestamp();
            double micros = Stopwatch.GetElapsedTime(start, end).TotalMicroseconds;
            return (MoleculeStatus.Rejected, Array.Empty<string>(), ex.Message, micros);
        }
    }

    // ------------------------------------------------------------------
    // Mode 1: name
    // ------------------------------------------------------------------

    private static void RunNameMode(Options options)
    {
        using var reader = new StreamReader(options.In!);
        using var writer = new StreamWriter(options.Out!);

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (string.IsNullOrWhiteSpace(line)) continue;

            var parsed = ParseLine(line);
            JsonObject result = parsed.Passthrough;

            if (parsed.PreValidationStatus.HasValue)
            {
                result["status"] = "invalid";
                result["error"] = parsed.PreValidationError;
                result["names"] = new JsonArray();
                result["micros"] = 0.0;
            }
            else
            {
                var (status, names, error, micros) = NameMolecule(parsed.Atoms!, parsed.Spec!.Value);
                result["status"] = status == MoleculeStatus.Named ? "named" : "rejected";
                result["names"] = new JsonArray(names.Select(n => (JsonNode)JsonValue.Create(n)).ToArray());
                result["micros"] = micros;
                if (error != null) result["error"] = error;
            }

            writer.WriteLine(result.ToJsonString());
        }
    }

    // ------------------------------------------------------------------
    // Mode 2: bench
    // ------------------------------------------------------------------

    private sealed class BenchMolecule
    {
        public required List<ElementGraph.AtomInput> Atoms;
        public required IUPAC.namingSpec Spec;
    }

    private static void RunBenchMode(Options options)
    {
        // Pre-deserialise everything first so JSON parsing is excluded from the timed region.
        var molecules = new List<BenchMolecule>();
        int invalidCount = 0;

        foreach (var line in File.ReadLines(options.In!))
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            var parsed = ParseLine(line);
            if (parsed.PreValidationStatus.HasValue)
            {
                invalidCount++;
                continue;
            }
            molecules.Add(new BenchMolecule { Atoms = parsed.Atoms!, Spec = parsed.Spec!.Value });
        }

        if (molecules.Count == 0)
        {
            throw new RunnerUsageException("No valid molecules found in input.");
        }

        // Warm-up passes (JIT warm-up), discarded.
        for (int w = 0; w < options.Warmup; w++)
        {
            foreach (var m in molecules)
            {
                NameMolecule(m.Atoms, m.Spec);
            }
        }

        // Timed passes. Single-threaded, run sequentially for honest per-molecule latency.
        var allLatenciesMicros = new List<double>(molecules.Count * options.Reps);
        var perRepSeconds = new double[options.Reps];
        int namedCount = 0;
        int rejectedCount = 0;

        var overallSw = Stopwatch.StartNew();
        for (int r = 0; r < options.Reps; r++)
        {
            long repStart = Stopwatch.GetTimestamp();
            foreach (var m in molecules)
            {
                var (status, _, _, micros) = NameMolecule(m.Atoms, m.Spec);
                allLatenciesMicros.Add(micros);
                if (r == options.Reps - 1)
                {
                    // Only tally named/rejected once (final rep) to avoid double counting;
                    // engine behaviour is deterministic across reps for the same input.
                    if (status == MoleculeStatus.Named) namedCount++;
                    else rejectedCount++;
                }
            }
            long repEnd = Stopwatch.GetTimestamp();
            perRepSeconds[r] = Stopwatch.GetElapsedTime(repStart, repEnd).TotalSeconds;
        }
        overallSw.Stop();

        double totalWallSeconds = overallSw.Elapsed.TotalSeconds;

        var throughputsPerRep = perRepSeconds.Select(s => molecules.Count / s).ToArray();
        double throughputMean = throughputsPerRep.Average();
        double throughputStdDev = StdDev(throughputsPerRep, throughputMean);

        allLatenciesMicros.Sort();
        double latMean = allLatenciesMicros.Average();
        double p50 = Percentile(allLatenciesMicros, 0.50);
        double p90 = Percentile(allLatenciesMicros, 0.90);
        double p99 = Percentile(allLatenciesMicros, 0.99);
        double p999 = Percentile(allLatenciesMicros, 0.999);
        double max = allLatenciesMicros[^1];

        var resultObj = new JsonObject
        {
            ["moleculeCount"] = molecules.Count,
            ["invalidInputCount"] = invalidCount,
            ["reps"] = options.Reps,
            ["warmup"] = options.Warmup,
            ["totalWallSeconds"] = totalWallSeconds,
            ["throughput"] = new JsonObject
            {
                ["moleculesPerSecMean"] = throughputMean,
                ["moleculesPerSecStdDev"] = throughputStdDev,
                ["perRepMoleculesPerSec"] = new JsonArray(throughputsPerRep.Select(t => (JsonNode)JsonValue.Create(t)).ToArray())
            },
            ["latencyMicros"] = new JsonObject
            {
                ["mean"] = latMean,
                ["p50"] = p50,
                ["p90"] = p90,
                ["p99"] = p99,
                ["p999"] = p999,
                ["max"] = max
            },
            ["outcomes"] = new JsonObject
            {
                ["named"] = namedCount,
                ["rejected"] = rejectedCount
            },
            ["machine"] = new JsonObject
            {
                ["osDescription"] = RuntimeInformation.OSDescription,
                ["processArchitecture"] = RuntimeInformation.ProcessArchitecture.ToString(),
                ["processorCount"] = Environment.ProcessorCount,
                ["frameworkDescription"] = RuntimeInformation.FrameworkDescription
            }
        };

        File.WriteAllText(options.Out!, resultObj.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));

        Console.WriteLine($"Benchmarked {molecules.Count} molecules x {options.Reps} reps " +
            $"({invalidCount} invalid inputs skipped). " +
            $"Throughput: {throughputMean:F1} mol/s (stddev {throughputStdDev:F1}). " +
            $"Latency p50={p50:F1}us p99={p99:F1}us max={max:F1}us.");
    }

    private static double StdDev(double[] values, double mean)
    {
        if (values.Length <= 1) return 0.0;
        double sumSq = values.Sum(v => (v - mean) * (v - mean));
        return Math.Sqrt(sumSq / (values.Length - 1));
    }

    private static double Percentile(List<double> sortedValues, double p)
    {
        if (sortedValues.Count == 0) return 0.0;
        double rank = p * (sortedValues.Count - 1);
        int lo = (int)Math.Floor(rank);
        int hi = (int)Math.Ceiling(rank);
        if (lo == hi) return sortedValues[lo];
        double frac = rank - lo;
        return sortedValues[lo] + (sortedValues[hi] - sortedValues[lo]) * frac;
    }
}
