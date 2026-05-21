# pprof and analyze_profile.py wrapper reference

This reference is loaded on demand by the `profile-and-optimize` skill
when reading profile output. The project's `scripts/analyze_profile.py`
wrapper is the right entry point for most analysis work — it knows about
Lancet2's layer structure and applies project-specific grouping. Reach
for raw pprof when the wrapper doesn't expose what you need.


Lancet2's profile analysis goes through `scripts/analyze_profile.py`, which wraps pprof with project-specific module/component classification, baseline tracking, and HTML rendering. For most analysis work the wrapper is the right entry point — it handles binary auto-detection, applies the project's grouping conventions, and produces structured rich-text reports. Reach for raw pprof when you need to do something the wrapper doesn't expose, or when you want to drive pprof's interactive shell for ad-hoc exploration.

The wrapper has its own help text:

```bash
pixi run -e profiling analyze-profile <profile.bin> -- --help
```

For pprof's full reference, see [github.com/google/pprof/blob/main/doc/README.md](https://github.com/google/pprof/blob/main/doc/README.md).

### Wrapper views

`analyze_profile.py` exposes named views that map onto pprof reports plus the project's classification layer:

| View | What it shows |
|:---|:---|
| `overview` | Total samples, sampling rate, top-level component breakdown |
| `top` | Flat and cumulative hottest functions (pprof's `-top`, with project demangling) |
| `modules` | Grouping by Lancet2 layer and external dependency |
| `components` | High-level component attribution (cbdg vs caller vs hts vs ...) |
| `tree` | Caller/callee tree, optionally focused with `--focus` |
| `hotpaths` | Ranked end-to-end call chains |
| `lines` | Source-line attribution (NOT included in `--view all`) |
| `all` | Everything except `lines` |

Useful invocations:

```bash
pixi run -e profiling analyze-profile profile.bin -- --view top --top 50
pixi run -e profiling analyze-profile profile.bin -- --view tree --focus 'lancet::caller'
pixi run -e profiling analyze-profile profile.bin -- --view hotpaths
pixi run -e profiling analyze-profile profile.bin -- --html /tmp/profile.html
pixi run -e profiling analyze-profile profile.bin -- --list 'BuildGraph'
pixi run -e profiling analyze-profile profile.bin -- --save-summary baseline
pixi run -e profiling analyze-profile -- --diff-tag baseline candidate1
pixi run -e profiling analyze-profile --history
```

The `--save-summary <name>` / `--diff-tag <a> <b>` flow is the project's canonical optimization workflow because it survives binary rebuilds. Plain `--diff-base` works only when both profiles came from the same binary.

### Raw pprof reference

For when the wrapper isn't enough. The basic shape is `pprof <format> [options] <binary> <profile>`:

```bash
# Locate pprof (the pixi profiling env installs it).
pixi run -e profiling ensure-pprof
PPROF=$(pixi run -e profiling which pprof)

# Interactive shell (no format flag).
$PPROF cmake-build-relwithdebinfo/Lancet2 Lancet.cpu_profile.<TS>.bin

# Web UI (browser-based, flame graph + graph + source).
$PPROF -http=localhost:8080 cmake-build-relwithdebinfo/Lancet2 Lancet.cpu_profile.<TS>.bin
```

Report formats (pick one):

| Flag | What it produces |
|:---|:---|
| `-text` / `-top` | One line per function with flat/cum samples |
| `-tree` | Callgraph in text, each function with its predecessors and successors |
| `-peek=<re>` | Full callers/callees of functions matching regex (no trimming) |
| `-traces` | Every sample as a stack trace, one location per line |
| `-list=<re>` | Source listing of functions matching regex, line-annotated |
| `-disasm=<re>` | Disassembly listing of functions matching regex |
| `-svg`, `-pdf`, `-png`, `-gif`, `-dot` | Graphical callgraph in the requested format |
| `-web` | Generates SVG and opens in a browser |
| `-weblist=<re>` | Combined source + assembly listing, browser-rendered |
| `-callgrind` | Callgrind format, for kcachegrind |

Common filtering options (combinable across most formats):

| Option | Effect |
|:---|:---|
| `-flat` (default) / `-cum` | Sort by flat or cumulative samples |
| `-functions` (default) / `-files` / `-lines` / `-addresses` | Granularity of report nodes |
| `-noinlines` | Attribute inlined functions to their out-of-line caller |
| `-nodecount=<n>` | Cap the report at n entries (default 80) |
| `-nodefraction=<f>` | Drop nodes below this fraction of total (default 0.005) |
| `-edgefraction=<f>` | Drop edges below this fraction of total (default 0.001) |
| `-focus=<re>` | Keep only call paths containing a node matching regex |
| `-ignore=<re>` | Drop call paths containing a node matching regex |
| `-show=<re>` / `-hide=<re>` | Show/hide specific entries |
| `-show_from=<re>` | Drop entries above the first match |

Comparison options:

| Option | Effect |
|:---|:---|
| `-diff_base=<profile>` | Subtract this profile from the source; percentages relative to base. Only valid for compatible same-binary profiles. |
| `-base=<profile>` | Subtract a cumulative base; percentages relative to the difference. |
| `-normalize` | Scale source profile total to match base before subtracting. Useful when comparing runs of different durations. |

Symbolization (rarely needs override; defaults work for Lancet2):

| Option | Effect |
|:---|:---|
| `-symbolize=local` | Symbolize from local binaries only |
| `-symbolize=demangle=full` | Demangle C++ names without parameter trimming |
| `-symbolize=demangle=templates` | Demangle, drop function params but keep template params |
| `-symbolize=none` | Disable symbolization (rare; for raw address profiles) |

### Interactive-shell commands

When started without a format flag, pprof drops into an interactive shell. The most useful commands match the report formats above (`top`, `tree`, `peek`, `list`, `disasm`, `traces`, `web`). Options can be set with `option=value`:

```
(pprof) top 30
(pprof) focus=lancet::caller
(pprof) cum
(pprof) top 30
(pprof) list BuildGraph
(pprof) peek BuildGraph
(pprof) web
```

Type `help` at the prompt for the full command list.

### When to use raw pprof vs the wrapper

The wrapper covers the common cases — top, tree, hotpaths, modules, components, source listing, save-summary, diff-tag, HTML report. Reach for raw pprof when you need:

- `pprof -http=localhost:8080` for the interactive web UI with flame graphs (the wrapper produces static HTML; the web UI lets you click through).
- `pprof -peek=<re>` for "show me everything around this one function" (no trimming, unlike `-tree`).
- `pprof -callgrind` for kcachegrind-based exploration.
- Custom `-focus` / `-ignore` / `-show_from` combinations the wrapper doesn't expose.
- `pprof -tagfocus`, `-tagignore`, `-tagroot`, `-tagleaf` — sample tag filtering for profiles that carry tags (gperftools profiles don't, but other profile sources do).

Otherwise, prefer the wrapper. Its module/component classification encodes Lancet2's layer structure, which raw pprof cannot know about.
