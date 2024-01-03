import numpy as np

rng = np.random.default_rng()

def fake_peak_trace(halfsize: int) -> np.ndarray:
    """
    Generates a random fake peak, with measurements increasing in value up to the center point and decreasing after.

    This currently operates with sorted exponential vectors, because that hack generates a visually pleasing peak. Should probably do something actually meaningful.
    
    TODO: There is obviously shared and repeated functionality between the peak methods; make things more modular
    """
    first_half = np.sort(rng.exponential(size=halfsize, scale=0.025))
    second_half = np.flip(np.sort(rng.exponential(size=halfsize, scale=0.025)))
    return np.concatenate([first_half, second_half])

def expspace_zero_one(num: int,
                      a: float) -> np.ndarray:
    """
    Return numbers spaced over the interval 0 to 1 along an exponential curve.

    Calculated as y = (a^x - 1) / (a - 1), a > 1.

    a controls the depth of the curve; higher values result in a longer wait before going to 1.
    """
    if a <= 1:
        raise ValueError('Value of a must be > 1.')
    return (np.power(a, np.linspace(start=0, stop=1, num=num)) - 1) / (a - 1)

def fake_read_mod_calls(halfsize: int,
                        read_type: str) -> np.ndarray:
    """
    Generates a read of the given size with modifications; returns 0 where there is no mod, 1 where there is a mod.

    TODO: More realistic read varieties
    """
    # Set the vector of p-vals for the bernoulli distribution pulls based on the requested read type
    match read_type:
        case 'peak':
            # higher chance of mod at center of read
            p_vec = expspace_zero_one(num=halfsize, a=15)
        case 'uniform':
            # uniform low chance of mod across entire read
            p_vec = [0.05] * halfsize
        case 'inverse_peak':
            # higher chance of mod at edges of read
            p_vec = np.flip(expspace_zero_one(num=halfsize, a=15))
        case _:
            ValueError(f'Unknown read type {read_type}')
    first_half = [np.random.binomial(n=1, p=x) for x in p_vec]
    second_half = np.flip([np.random.binomial(n=1, p=x) for x in p_vec])
    return np.concatenate([first_half, second_half])

def fake_read_mod_positions(halfsize: int,
                            read_type: str) -> np.ndarray:
    """
    Generates a read of the given size with modifications; returns positions where there is a modification.

    Positions are relative to the center of the read.

    See fake_read_mod_calls for details.
    """
    return np.flatnonzero(fake_read_mod_calls(halfsize=halfsize, read_type=read_type)) - halfsize