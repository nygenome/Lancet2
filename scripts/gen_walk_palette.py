#!/usr/bin/env python3
"""Generate maximally-distinct colors for walk edge visualization in Graphviz DOT.

Produces a pre-computed C++ constexpr lookup table of 64 hex color strings that
are maximally distinguishable to the human eye.

Algorithm (inspired by iwanthue — https://medialab.github.io/iwanthue/theory/):
  1. SAMPLE in sRGB [0,1]^3 — guaranteed in-gamut, no clipping warnings.
  2. CONVERT to CIE L*a*b* — the perceptually uniform color space where
     Euclidean distance approximates perceived color difference (Delta-E).
  3. FILTER by perceptual constraints in L*a*b*:
       L* ∈ [30, 75]  — avoid near-black / near-white (poor graph visibility).
       Chroma ∈ [35, 110] — avoid grays (indistinguishable) and neon (jarring).
  4. CLUSTER via k-means in L*a*b* — finds N maximally-separated centroids.
  5. REORDER via farthest-first traversal — ensures any prefix [0..K) has
     maximum pairwise separation. This is critical because walk indices are
     assigned sequentially (walk 1, walk 2, ...) and most components have
     only a few walks. Adjacent indices should look maximally different.
  6. CONVERT centroids back to sRGB hex (#RRGGBB) for Graphviz DOT.

Why sample RGB → filter in L*a*b* (not sample L*a*b* directly)?
  The sRGB gamut is a twisted irregular volume inside L*a*b*. Sampling random
  L*a*b* points produces many out-of-gamut colors (negative XYZ values) that
  must be clipped or rejected. Sampling RGB first guarantees every point is a
  valid displayable color. The L*a*b* conversion is then used only for the
  perceptual filter and the k-means distance metric.

Usage:
  pixi run -e walk-palette python scripts/gen_walk_palette.py
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.spatial.distance import pdist, squareform
from skimage.color import lab2rgb, rgb2lab
from sklearn.cluster import MiniBatchKMeans

# ── L*a*b* Perceptual Constraints ───────────────────────────────────────────
# L* ∈ [30, 75]: skip near-black and near-white for graph background contrast.
# Chroma (sqrt(a*² + b*²)) ∈ [35, 110]: skip desaturated grays and extreme neon.
L_STAR_MIN, L_STAR_MAX = 30, 75
CHROMA_MIN, CHROMA_MAX = 35, 110

NUM_COLORS = 64
NUM_RGB_SAMPLES = 500_000  # dense sRGB sampling for good L*a*b* coverage
RANDOM_SEED = 42


def sample_constrained_lab(num_samples: int, seed: int) -> np.ndarray:
    """Sample sRGB colors uniformly, convert to L*a*b*, apply constraints.

    Returns an (M, 3) array of L*a*b* points that pass the perceptual filter,
    where M <= num_samples (some are rejected by the chroma/lightness filter).
    """
    rng = np.random.default_rng(seed)
    rgb = rng.uniform(0.0, 1.0, size=(num_samples, 3))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lab = rgb2lab(rgb.reshape(-1, 1, 3)).reshape(-1, 3)

    l_star = lab[:, 0]
    chroma = np.sqrt(lab[:, 1] ** 2 + lab[:, 2] ** 2)
    mask = (
        (l_star >= L_STAR_MIN) & (l_star <= L_STAR_MAX)
        & (chroma >= CHROMA_MIN) & (chroma <= CHROMA_MAX)
    )

    filtered = lab[mask]
    print(f"// Sampling: {num_samples} RGB -> {len(filtered)} pass L*a*b* filter "
          f"({100 * len(filtered) / num_samples:.1f}% yield)")
    return filtered


def cluster_palette(lab_samples: np.ndarray, n_colors: int, seed: int) -> np.ndarray:
    """Run MiniBatchKMeans in L*a*b* space to find n_colors cluster centers.

    MiniBatchKMeans is ~10x faster than KMeans for large sample counts with
    negligible quality loss on well-separated clusters like color data.
    """
    kmeans = MiniBatchKMeans(
        n_clusters=n_colors,
        random_state=seed,
        n_init=10,
        max_iter=300,
        batch_size=min(10_000, len(lab_samples)),
    )
    kmeans.fit(lab_samples)
    return kmeans.cluster_centers_


def farthest_first_order(centers: np.ndarray) -> np.ndarray:
    """Reorder centers so each successive color is maximally distant from all
    previously selected colors. This ensures any prefix [0..K) has the best
    possible minimum pairwise separation.

    Greedy algorithm (O(N^2)):
      1. Start with the color farthest from the centroid of all colors.
      2. At each step, pick the unselected color whose minimum distance to
         any already-selected color is largest (maximize the bottleneck).
    """
    n = len(centers)
    dist_matrix = squareform(pdist(centers))

    # Start with color farthest from the centroid
    centroid = centers.mean(axis=0)
    dists_to_centroid = np.linalg.norm(centers - centroid, axis=1)
    first = int(np.argmax(dists_to_centroid))

    selected = [first]
    remaining = set(range(n)) - {first}

    # min_dist_to_selected[i] = min distance from i to any selected color
    min_dist_to_selected = dist_matrix[first].copy()

    for _ in range(n - 1):
        # Pick the remaining color with the largest min distance to selected set
        best = max(remaining, key=lambda i: min_dist_to_selected[i])
        selected.append(best)
        remaining.remove(best)

        # Update min distances
        np.minimum(min_dist_to_selected, dist_matrix[best], out=min_dist_to_selected)

    return centers[selected]


def lab_to_hex(lab_points: np.ndarray) -> list[str]:
    """Convert L*a*b* array (N, 3) to list of #RRGGBB hex strings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rgb = lab2rgb(lab_points.reshape(-1, 1, 3)).reshape(-1, 3)
    rgb = np.clip(rgb, 0.0, 1.0)
    return [f"#{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}"
            for r, g, b in rgb]


def compute_quality_metrics(centers: np.ndarray, hex_colors: list[str]) -> None:
    """Print pairwise Delta-E statistics for the full palette and prefix subsets."""
    all_dists = pdist(centers)
    full_sq = squareform(all_dists)
    np.fill_diagonal(full_sq, np.inf)

    min_de = np.min(full_sq)
    min_i, min_j = np.unravel_index(np.argmin(full_sq), full_sq.shape)

    print(f"//")
    print(f"// Delta-E quality (Euclidean L*a*b* distance):")
    print(f"//   dE > 10 is easily distinguishable, dE > 20 is very distinct.")
    print(f"//   Full palette: min dE = {min_de:.1f} "
          f"([{min_i}]={hex_colors[min_i]} vs [{min_j}]={hex_colors[min_j]}), "
          f"avg dE = {np.mean(all_dists):.1f}")

    # Fine-grained prefixes for small N (most common case), coarser for large N
    prefixes = list(range(2, 21)) + [24, 32, 48, 64]
    for n in prefixes:
        if n > len(centers):
            break
        sub = centers[:n]
        sub_sq = squareform(pdist(sub))
        np.fill_diagonal(sub_sq, np.inf)
        sub_min = np.min(sub_sq)
        si, sj = np.unravel_index(np.argmin(sub_sq), sub_sq.shape)
        print(f"//   First {n:2d}: min dE = {sub_min:5.1f} "
              f"([{si}]={hex_colors[si]} vs [{sj}]={hex_colors[sj]})")


def emit_cpp_array(hex_colors: list[str], centers: np.ndarray) -> None:
    """Print C++ constexpr array definition with per-element L*a*b* metadata."""
    print(f"// {len(hex_colors)} maximally-distinct colors for walk edge visualization.")
    print(f"// Generated by scripts/gen_walk_palette.py via k-means in CIE L*a*b*.")
    print(f"// Constraints: L* in [{L_STAR_MIN}, {L_STAR_MAX}], "
          f"Chroma in [{CHROMA_MIN}, {CHROMA_MAX}].")
    compute_quality_metrics(centers, hex_colors)
    print(f"//")
    print(f"// Ordered by farthest-first traversal: any prefix [0..K) has")
    print(f"// maximum pairwise separation, critical for small walk counts.")

    n = len(hex_colors)
    print(f"static constexpr std::array<std::string_view, {n}> WALK_COLORS = {{")
    for i, (hx, lab) in enumerate(zip(hex_colors, centers)):
        comma = "," if i < n - 1 else " "
        hue = np.degrees(np.arctan2(lab[2], lab[1])) % 360
        chroma = np.sqrt(lab[1] ** 2 + lab[2] ** 2)
        print(f'    "{hx}"sv{comma}  '
              f"// [{i:2d}] L*={lab[0]:4.1f} C={chroma:4.1f} H={hue:5.1f}")
    print("};")


def main() -> None:
    lab_samples = sample_constrained_lab(NUM_RGB_SAMPLES, RANDOM_SEED)
    centers = cluster_palette(lab_samples, NUM_COLORS, RANDOM_SEED)
    centers = farthest_first_order(centers)
    hex_colors = lab_to_hex(centers)
    emit_cpp_array(hex_colors, centers)


if __name__ == "__main__":
    main()
