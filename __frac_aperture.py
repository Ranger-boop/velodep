import ctypes
import numpy as np
from .__load_dll import *


# Define the argtypes and restypes of the C apis.
dll.new_FracAperture.restype = ctypes.POINTER(ctypes.c_void_p)
dll.new_FracAperture.argtypes = (
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)

dll.delete_FracAperture.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.solve_FracAperture.restype = None
dll.solve_FracAperture.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_len_b_FracAperture.restype = ctypes.c_size_t
dll.get_len_b_FracAperture.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_b_mod_FracAperture.restype = ctypes.POINTER(ctypes.c_double)
dll.get_b_mod_FracAperture.argtypes = (ctypes.POINTER(ctypes.c_void_p),)


class FracAperture:
    """
    A solver to compute fracture aperture based on the displacement- and
    velocity-dependent aperture model.

    Parameters
    ----------
    dil_fact : float
        Dilation factor.
    D_c : float
        Characteristic slip distance.
    dil_ang : float
        Dilation angle.
    b_0 : float
        Initial aperture at the onset of reactivation (the beginning of the
        first velocity step).
    u_end : ndarray
        Slip displacement after the reactivation of fracture at the end of
        each velocity step.
    u_ini : ndarray
        Slip displacement after the reactivation of fracture at the initial
        time of each velocity step.
    v : ndarray
        Slip velocity at each velocity step.
    cp : bool
        A flag to specify if it is necessary to copy the transfered arrays.
        The default value is True for memory safety. Turn it to False only if
        the transfered arrays are self-holded np.ndarray objects rather than
        the slices of any other objects.
    """

    def __init__(
        self,
        dil_fact: np.float64,
        D_c: np.float64,
        dil_ang: np.float64,
        b_0: np.float64,
        u_end: np.ndarray,
        u_ini: np.ndarray,
        v: np.ndarray,
        cp=True,
    ) -> None:
        if cp:
            self.u_end = u_end.copy()
            self.u_ini = u_ini.copy()
            self.v = v.copy()
        else:
            self.u_end = u_end
            self.u_ini = u_ini
            self.v = v
        ptr_u_end = self.u_end.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_u_end = self.u_end.size
        ptr_u_ini = self.u_ini.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_u_ini = self.u_ini.size
        ptr_v = self.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_v = self.v.size
        self.obj = dll.new_FracAperture(
            dil_fact,
            D_c,
            dil_ang,
            b_0,
            ptr_u_end,
            size_u_end,
            ptr_u_ini,
            size_u_ini,
            ptr_v,
            size_v,
        )

    def __del__(self) -> None:
        dll.delete_FracAperture(self.obj)

    def solve(self) -> None:
        dll.solve_FracAperture(self.obj)

    @property
    def __len_b(self) -> ctypes.c_size_t:
        return dll.get_len_b_FracAperture(self.obj)

    @property
    def b_mod(self) -> np.ndarray:
        ptr = dll.get_b_mod_FracAperture(self.obj)
        return np.ctypeslib.as_array(ptr, shape=(self.__len_b,))
