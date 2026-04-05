# Graph Complexity (`GRAPH_CX`)

Graph-level metrics describe the structural complexity of the colored de
Bruijn graph component surrounding a variant. Requires
`--enable-graph-complexity-features` to enable.

---

## Field Layout

**INFO tag**: `GRAPH_CX` (3 comma-separated values)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 0 | GEI | Float | [0, ∞) | Graph Entanglement Index — log-squashed topology chaos metric |
| 1 | TipToPathCovRatio | Float | [0, ∞) | Mean dead-end coverage / mean unitig coverage |
| 2 | MaxSingleDirDegree | Integer | [0, ∞) | Maximum outgoing edges in any single sign direction |

---

## Graph Entanglement Index (GEI)

The GEI compresses four collinear topology metrics into a single continuous
scalar. Before defining the formula, here is what each input measures:

### GEI Input Metrics

#### Cyclomatic Complexity (CC)

- **Intuition**: How many independent loops exist in the graph. A linear
  chain of k-mers has CC=0. A clean heterozygous SNP creates exactly one
  bubble (CC=1). A tandem repeat region creates many nested loops.
- **Math**: `CC = E − V + 1` for a single connected component, where E is
  the number of edges and V is the number of nodes. This is the graph-theoretic
  circuit rank — the minimum number of edges you must remove to make the
  graph a tree.
- **Biology**: CC=0–1 means the assembler found one clean path divergence
  (a normal variant). CC=5–10 means overlapping repeat units are creating
  multiple competing paths. CC>50 typically means a tandem repeat region
  (STR/VNTR) or segmental duplication where many near-identical k-mers
  create a hairball of alternative paths.

#### Branch Points (BP)

- **Intuition**: How many crossroads the assembler encounters. Each branch
  point is a node where the graph splits into ≥2 outgoing edges in at least
  one direction — a "fork in the road" for path enumeration.
- **Math**: Count of nodes with `max(out_degree_fwd, out_degree_rev) ≥ 2`.
  In a clean variant bubble, exactly 2 branch points exist (the bubble's
  entry and exit nodes). In pathological graphs, nearly every node branches.
- **Biology**: A clean heterozygous variant produces BP=2 (one fork, one
  merge). BP≥5 in a ~1000bp window means multiple overlapping variant
  signals. BP≥50 indicates the region is dominated by repetitive sequence
  where the k-mer size is too short to resolve unique paths — the graph
  is "shattered" into many short fragments connected by ambiguous k-mers.

#### Coverage Coefficient of Variation (CoverageCv)

- **Intuition**: How uniformly reads are distributed across the graph. If
  every node has roughly the same number of supporting reads, CovCV is low.
  If some nodes have 500× coverage while neighbors have 5×, CovCV is high.
- **Math**: `CovCV = σ / μ` — the standard deviation of total read support
  across all nodes divided by the mean. Dimensionless ratio where CovCV=0
  means perfectly uniform coverage and CovCV>1 means the standard deviation
  exceeds the mean.
- **Biology**: True biological variants produce near-uniform coverage
  (CovCV<0.5) because both haplotypes are sequenced at similar depth.
  Collapsed repeats produce extreme coverage spikes (CovCV>1.5) because
  reads from 2+ genomic loci pile onto the same graph nodes. This is the
  key distinguisher between "complex but real" (high CC, low CovCV) and
  "collapsed artifact" (high CC, high CovCV).

#### Unitig Ratio

- **Intuition**: What fraction of the graph is "boring" linear chain. A
  unitig node has exactly one incoming and one outgoing edge — there is no
  ambiguity about what comes before or after it. High unitig ratio means
  most of the graph is resolved; low means it's fracturing everywhere.
- **Math**: `UnitigRatio = (nodes with in_degree=1 and out_degree=1) / V`.
  Ranges from 0.0 (every node branches) to 1.0 (perfect linear chain;
  impossible in practice since at least one source and one sink exist).
- **Biology**: Clean variant regions have UnitigRatio>0.90 — the vast
  majority of k-mers are unambiguous, with just a few branch points at
  the variant boundaries. STR regions drop to 0.60–0.80. Paralogous
  collapses can reach 0.10–0.30, meaning 70–90% of all positions in the
  graph are ambiguous forks.

### Formula

```
GEI = log₁₀(1 + (CC × BP × CoverageCv) / (UnitigRatio + ε))
```

**Why compress?** These metrics are mathematically interlocked by graph
theory identities. Increasing BP necessarily decreases UnitigRatio and
increases CC. Feeding all four to an additive ML model causes "vote
splitting" — the model distributes signal weight across correlated features,
producing noisy, low-confidence decision curves.

**The multiplicative AND logic**: The multiplication acts as a soft AND gate.
All three numerator terms must be elevated for the product to spike. A graph
with CC=100 but BP=1 and CovCV=0.1 yields a low GEI — it's complex on
one axis but not entangled.

**Component roles**:

- **Numerator** (CC × BP × CoverageCv): Drives the chaos signal.
  CoverageCv distinguishes true polymorphism (uniform coverage, low CV)
  from collapsed repeats (reads stacking on hub nodes, high CV).
- **Denominator** (UnitigRatio + ε): High UnitigRatio (~0.95) → no
  amplification. Low UnitigRatio (~0.10) → 10× multiplier, penalizing
  shattered graphs.
- **Log₁₀ squash**: Compresses raw values (which can be arbitrarily large)
  to a practical scale.

### Biological Interpretation

| GEI Range | Meaning |
|:----------|:--------|
| < 0.5 | **Pristine** — true biological SNP/INDEL in unique sequence. Clean bubble. |
| 0.5–2.0 | **STR / Microsatellite** — k-mer struggles to bridge exact repeats, creating local loops, but flanking anchors hold. |
| 2.0–3.0 | **Complex region** — multiple overlapping events, possible SV breakpoints. |
| > 3.0 | **Paralogous collapse** — segmental duplication or centromeric satellite. Reads from multiple genomic loci collapse into the same graph. Any variant called here is almost certainly an artifact. |

---

## TipToPathCovRatio

Ratio of mean dead-end (tip) node coverage to mean linear (unitig) node
coverage. Measures assembly tearing and structural variant breakpoints,
independent of loop topology.

| Range | Meaning |
|:------|:--------|
| ≈ 0 | Clean — tips have negligible coverage |
| 0.5–2.0 | Moderate — some assembly fragmentation |
| > 2.0 | High — assembly is tearing apart, likely SV or mapping artifact |

---

## MaxSingleDirDegree

Maximum outgoing edges in any single sign direction across all nodes in
the component. Detects hub k-mers (e.g., `AAAAAAA`) that act as black
holes, pulling in many unrelated graph paths.

| Range | Meaning |
|:------|:--------|
| ≤ 3 | Normal |
| 4–5 | Elevated — possible low-complexity region |
| > 5 | Hub k-mer — likely homopolymer or simple repeat causing BFS blowup |
