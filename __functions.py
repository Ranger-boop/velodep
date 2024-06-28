import numpy as np


def check_0(vlc: np.ndarray) -> None:
    """
    Check how many zero(s) are in the vlocities.

    Parameters
    ----------
    vlc : ndarray
        The ndarray of velocities need to be checked.
    output : text.
        The text displaying the number of zero(s).
    """
    num_0 = 0
    for i in range(len(vlc)):
        if vlc[i] == 0:
            num_0 += 1
    print(f"{num_0} 0s are in the velocity list.")


def rmse(val_mod: np.ndarray, val_exp: np.ndarray) -> float:
    """
    Compute and return the value of Root Mean Square Error (RMSE)
    between the calculated and measured values.

    Parameters
    ----------
    val_mod : ndarray
        The values calculated by model.
    disp_exp : ndarray
        The values measured in experiment.
    output : float
        The computed RMSE of the two sets of values.
    """
    if np.isnan(val_mod).any():
        dest = np.inf
    else:
        dest = np.sqrt(np.sum((val_mod - val_exp) ** 2) / len(val_mod))
    return dest
