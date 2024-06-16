import os
import ctypes
import numpy as np

# Load the shared library
current_dir = os.path.dirname(os.path.abspath(__file__))
# dll_path = os.path.join(current_dir, "build/bin/Debug/rsf.dll")
dll_path = os.path.join(current_dir, "build/bin/Release/rsf.dll")
dll = ctypes.cdll.LoadLibrary(dll_path)

# Define the argtypes and restypes of the C apis.
dll.new_RSF_tau_fix.restype = ctypes.POINTER(ctypes.c_void_p)
dll.new_RSF_tau_fix.argtypes = (
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
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

dll.delete_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.solve_RSF.restype = None
dll.solve_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_len_RSF_sol.restype = ctypes.c_size_t
dll.get_len_RSF_sol.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_disp_RSF.restype = ctypes.POINTER(ctypes.c_double)
dll.get_disp_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_theta_RSF.restype = ctypes.POINTER(ctypes.c_double)
dll.get_theta_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_tau_RSF.restype = ctypes.POINTER(ctypes.c_double)
dll.get_tau_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_vel_RSF.restype = ctypes.POINTER(ctypes.c_double)
dll.get_vel_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

dll.get_mu_RSF.restype = ctypes.POINTER(ctypes.c_double)
dll.get_mu_RSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)


class RSF_tau_fix:
    """
    Runge-Kutta solver for the equation group of the rate-and-state friction law
    with the shear stress on the fault fixed.

    The solver is based on the 4th order Runge-Kutta method.

    Parameters
    ----------
    a : float
        Scaling factor of the direct velocity effect on friction coefficient.
    b : float
        Scaling factor of the evolutional effect on friction coefficient.
    Dc : float
        Characteristic distance over which the contact population has been renewed.
    alpha : float
        Effect of changing effective normal stress on shear stress.
    mu_0 : float
        Friction coefficient at a reference slip velocity of v_0.
    v_0 : float
        Reference slip velocity corresponding to mu_0.
    disp_0: float
        Initial displacement at the slip velocity of v_0.
    sig_n : numpy.ndarray
        Normal stress on the fault.
    pp : numpy.ndarray
        Pore pressure on the fault.
    time: numpy.ndarray
        Time sequence of normal stress and pore pressure.
    cp : bool
        A flag for whether or not to copy the numpy.ndarray objects.
        The default value is True to avoid pointer misalignment.
        Set it as False only if the numpy.ndarray objects are not slices of
        any data but original numpy.ndarray objects.
    """

    def __init__(
        self,
        a: np.float64,
        b: np.float64,
        Dc: np.float64,
        alpha: np.float64,
        mu_0: np.float64,
        v_0: np.float64,
        disp_0: np.float64,
        sig_n: np.ndarray,
        pp: np.ndarray,
        time: np.ndarray,
        cp=True,
    ) -> None:
        if cp:
            self.sig_n = sig_n.copy()
            self.pp = pp.copy()
            self.time = time.copy()
        else:
            self.sig_n = sig_n
            self.pp = pp
            self.time = time
        ptr_sig_n = self.sig_n.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_sig_n = self.sig_n.size
        ptr_pp = self.pp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_pp = self.pp.size
        ptr_time = self.time.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_time = self.time.size
        self.obj = dll.new_RSF_tau_fix(
            a,
            b,
            Dc,
            alpha,
            mu_0,
            v_0,
            disp_0,
            ptr_sig_n,
            size_sig_n,
            ptr_pp,
            size_pp,
            ptr_time,
            size_time,
        )

    def __del__(self) -> None:
        dll.delete_RSF(self.obj)

    def solve_rk4(self) -> None:
        dll.solve_RSF(self.obj)

    @property
    def len_RSF_sol(self) -> ctypes.c_size_t:
        return dll.get_len_RSF_sol(self.obj)

    @property
    def disp(self) -> np.ndarray:
        ptr = dll.get_disp_RSF(self.obj)
        return np.ctypeslib.as_array(ptr, shape=(self.len_RSF_sol,))

    @property
    def theta(self) -> np.ndarray:
        ptr = dll.get_theta_RSF(self.obj)
        return np.ctypeslib.as_array(ptr, shape=(self.len_RSF_sol,))

    @property
    def tau(self) -> np.ndarray:
        ptr = dll.get_tau_RSF(self.obj)
        return np.ctypeslib.as_array(ptr, shape=(self.len_RSF_sol,))

    @property
    def vel(self) -> np.ndarray:
        ptr = dll.get_vel_RSF(self.obj)
        return np.ctypeslib.as_array(ptr, shape=(self.len_RSF_sol,))

    @property
    def mu(self) -> np.ndarray:
        ptr = dll.get_mu_RSF(self.obj)
        return np.ctypeslib.as_array(ptr, shape=(self.len_RSF_sol,))


# Define the argtypes and restypes of the C apis.
dll.new_RSF_tau_chg.restype = ctypes.POINTER(ctypes.c_void_p)
dll.new_RSF_tau_chg.argtypes = (
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
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
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)


class RSF_tau_chg(RSF_tau_fix):
    """
    Runge-Kutta solver for the equation group of the rate-and-state friction law
    with the shear stress on the fault changable.

    The solver is based on the 4th order Runge-Kutta method.

    Parameters
    ----------
    a : float
        Scaling factor of the direct velocity effect on friction coefficient.
    b : float
        Scaling factor of the evolutional effect on friction coefficient.
    Dc : float
        Characteristic distance over which the contact population has been renewed.
    alpha : float
        Effect of changing effective normal stress on shear stress.
    k_s : float
        Stiffness of medium around the fault.
    mu_0 : float
        Friction coefficient at a reference slip velocity of v_0.
    v_0 : float
        Reference slip velocity corresponding to mu_0.
    disp_0: float
        Initial displacement at the slip velocity of v_0.
    vlp : ndarray
        Load point velocity (velocity of the medium around the fault).
    sig_n : ndarray
        Normal stress on the fault.
    pp : ndarray
        Pore pressure on the fault.
    time: ndarray
        Time sequence of normal stress and pore pressure.
    """

    def __init__(
        self,
        a: np.float64,
        b: np.float64,
        Dc: np.float64,
        alpha: np.float64,
        k_s: np.float64,
        mu_0: np.float64,
        v_0: np.float64,
        disp_0: np.float64,
        vlp: np.ndarray,
        sig_n: np.ndarray,
        pp: np.ndarray,
        time: np.ndarray,
        cp=True,
    ) -> None:
        if cp:
            self.vlp = vlp.copy()
            self.sig_n = sig_n.copy()
            self.pp = pp.copy()
            self.time = time.copy()
        else:
            self.vlp = vlp
            self.sig_n = sig_n
            self.pp = pp
            self.time = time
        ptr_vlp = self.vlp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_vlp = self.vlp.size
        ptr_sig_n = self.sig_n.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_sig_n = self.sig_n.size
        ptr_pp = self.pp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_pp = self.pp.size
        ptr_time = self.time.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        size_time = self.time.size
        self.obj = dll.new_RSF_tau_chg(
            a,
            b,
            Dc,
            alpha,
            k_s,
            mu_0,
            v_0,
            disp_0,
            ptr_vlp,
            size_vlp,
            ptr_sig_n,
            size_sig_n,
            ptr_pp,
            size_pp,
            ptr_time,
            size_time,
        )


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
