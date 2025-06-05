import numpy as np

from . import utils


def expspace_prob(num: int, a: float, b: float = 1) -> np.ndarray:
    """
    Return probability values spaced over the interval 0 to b along an exponential curve.

    Calculated as y = ((a^x - 1) / (a - 1)) * b, a > 1, 0 < b <= 1

    Args:
        num: total length of space to return; same as num argument to np.linspace
        a: controls the depth of the curve; higher values result in a longer wait before going to 1. Must be >1.
        b: controls the max value of the curve. Must be between (0, 1].

    Returns:
        Array of probability values betweeen 0 and b, spaced along an exponential curve
    """
    if a <= 1:
        raise ValueError("Value of a must be > 1.")
    if b > 1:
        raise ValueError("Value of b must be between (0, 1].")
    return (np.power(a, np.linspace(start=0, stop=1, num=num)) - 1) / (a - 1) * b


def fake_read_mod_calls(halfsize: int, read_type: str, max_prob: float) -> np.ndarray:
    """
    Generates a read of the given size with modifications; returns 0 where there is no mod, 1 where there is a mod.

    TODO: More realistic read varieties
    TODO: Fewer magic numbers

    Args:
        halfsize: specifies length of output trace; final length will be 2*halfsize
        read_type: string name of desired read type; see match statement for available types
        max_prob: maximum probability of any single position being called as modified

    Returns:
        Array of 0s and 1s, patterned appropriately
    """
    # Set the vector of p-vals for the bernoulli distribution pulls based on the requested read type
    match read_type:
        case "peak":
            # higher chance of mod at center of read
            p_vec = expspace_prob(num=halfsize, a=15, b=max_prob)
        case "uniform":
            # uniform low chance of mod across entire read
            p_vec = [0.05] * halfsize
        case "inverse_peak":
            # higher chance of mod at edges of read
            p_vec = np.flip(expspace_prob(num=halfsize, a=15, b=max_prob))
        case _:
            ValueError(f"Unknown read type {read_type}")
    first_half = [utils.rng.binomial(n=1, p=x) for x in p_vec]
    second_half = np.flip([utils.rng.binomial(n=1, p=x) for x in p_vec])
    return np.concatenate([first_half, second_half])


def fake_read_mod_positions(
    halfsize: int, read_type: str, max_prob: float
) -> np.ndarray:
    """
    Generates a read of the given size with modifications; returns positions where there is a modification.
    Positions are relative to the center of the read.

    See fake_read_mod_calls for details.

    TODO: Should this be able to operate in a non-centered fashion?

    Returns:
        Vector of positions of modified bases
    """
    return (
        np.flatnonzero(
            fake_read_mod_calls(
                halfsize=halfsize, read_type=read_type, max_prob=max_prob
            )
        )
        - halfsize
    )


def fake_peak_enrichment_profile(
    halfsize: int, peak_height: float, n_reads: int = 100
) -> np.ndarray:
    """
    Generates a random fake peak, with measurements increasing in value up to the center point and decreasing after.
    Returns enrichment values (fraction modified bases) at each position.

    Args:
        halfsize: specifies length of output trace; final length will be 2*halfsize
        peak_height: max height of peak. Must be (0, 1].
        n_reads: number of reads to generate

    Returns:
        Array of enrichment values between 0 and 1, peaking at the middle
    """
    reads = [fake_read_mod_calls(halfsize, "peak", peak_height) for _ in range(n_reads)]
    modified_base_counts = np.sum(reads, axis=0)
    modified_fractions = np.divide(modified_base_counts, n_reads)
    return modified_fractions


def fake_peak_enrichment(
    halfsize: int, peak_height: float, n_reads: int = 100
) -> tuple[int, int]:
    """
    Generates total modification enrichment counts for a set of fake peak reads spanning some region.
    Returns enrichment values (fraction modified bases) summed across the entire region.

    Args:
        halfsize: specifies length of output trace; final length will be 2*halfsize
        peak_height: max height of peak. Must be (0, 1].
        n_reads: number of reads to generate

    Returns:
        tuple containing counts of (modified bases, total_bases)
    """
    reads = [fake_read_mod_calls(halfsize, "peak", peak_height) for _ in range(n_reads)]
    modified_bases = np.sum(reads)
    total_bases = np.sum(np.fromiter((len(read) for read in reads), dtype=int))
    return (modified_bases, total_bases)
