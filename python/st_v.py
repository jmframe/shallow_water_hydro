import matplotlib.pyplot as plt
import imageio
import yaml
import random
import xarray as xr
import numpy as np

class SAINT_VENANT:
    def __init__(self, gcs_key=None):
        """This is a dynamic hydraulics model to solve the water surface elevation over time.
        Args:
        """
        self.data_dir = None
        self.verbosity = 0
        self.vertical_boundary = None

    def initialize(self, config_file: str):
        """
        This function should perform all tasks that are to take place before entering the model’s time loop.
        Models should be refactored, if necessary, to read their inputs
        (which could include filenames for other input files) from a configuration file.
        Args:
            config_file (str): [description]. Defaults to None.
        """

        # -------------- Read in the configurations -------------------------#
        with open(config_file, "r") as fp:
            cfg = yaml.safe_load(fp)
        self._parse_config(cfg)

        # ========================================== #
        if self.data_dir is None:
            print("data dir is none, setting synthetic test case for initial head")
            self.set_synthetic_initial_head()
        else:
            print("Reading head and geometry files")
            self.read_initial_head_from_file()

        # ========================================== #
        self.tau = self.dt / self.dx / 2
        self.time = np.arange(0, self.tmax, self.dt)
        self.nt = self.time.shape[0]
        self.U = np.zeros((3, self.nx, self.ny))
        self.Us = np.zeros((3, self.nx, self.ny))
        self.Uss = np.zeros((3, self.nx, self.ny))
        self.F = np.zeros((3, self.nx, self.ny))
        self.E = np.zeros((3, self.nx, self.ny))
        self.lmda = np.zeros((3, self.nx, self.ny))
        self.S = np.zeros((3, self.nx, self.ny))
        self.v = np.zeros((self.nx, self.ny))
        self.u = np.zeros((self.nx, self.ny))
        self.x0 = np.round(self.nx / 2, 2)
        self.alt = None
        self.images = []

    # ========================================== #
    def predictor(self):
        """a "provisional" value of u at time level n+1
        (denoted by {\displaystyle u_{i}^{\overline {n+1}}}u_i^{\overline{n+1}}) is estimated as follows
        u_i^{\overline{n+1}} = u_i^n - a \frac{\Delta t}{\Delta x} \left( u_{i+1}^n - u_i^n \right)
        NOTE: There is a masking of the gradients if there is adjacent bare ground.
        """
        h = self.h
        dem = self.dem
        U = self.U
        Us = self.Us
        E = self.E
        F = self.F
        S = self.S
        tau = self.tau
        dt = self.dt
        alt = self.alt

        if alt == 0:
            E1 = np.s_[2:, 1:-1]
            E2 = np.s_[1:-1, 1:-1]
            F1 = np.s_[1:-1, 1:-1]
            F2 = np.s_[1:-1, :-2]
            Ediff = E[:, 2:, 1:-1] - E[:, 1:-1, 1:-1]
            Fdiff = F[:, 1:-1, 1:-1] - F[:, 1:-1, :-2]
        elif alt == 1:
            E1 = np.s_[1:-1, 1:-1]
            E2 = np.s_[:-2, 1:-1]
            F1 = np.s_[1:-1, 2:]
            F2 = np.s_[1:-1, 1:-1]
            Ediff = E[:, 1:-1, 1:-1] - E[:, :-2, 1:-1]
            Fdiff = F[:, 1:-1, 2:] - F[:, 1:-1, 1:-1]
        elif alt == 2:
            E1 = np.s_[2:, 1:-1]
            E2 = np.s_[1:-1, 1:-1]
            F1 = np.s_[1:-1, 2:]
            F2 = np.s_[1:-1, 1:-1]
            Ediff = E[:, 2:, 1:-1] - E[:, 1:-1, 1:-1]
            Fdiff = F[:, 1:-1, 2:] - F[:, 1:-1, 1:-1]
        elif alt == 3:
            E1 = np.s_[1:-1, 1:-1]
            E2 = np.s_[:-2, 1:-1]
            F1 = np.s_[1:-1, 1:-1]
            F2 = np.s_[1:-1, :-2]
            Ediff = E[:, 1:-1, 1:-1] - E[:, :-2, 1:-1]
            Fdiff = F[:, 1:-1, 1:-1] - F[:, 1:-1, :-2]

        Ediff = np.where(
            Ediff > 0,
            Ediff * ((h[E1] - dem[E1]) / (h[E1] + 1)),
            Ediff * ((h[E2] - dem[E2]) / (h[E2] + 1)),
        )
        Fdiff = np.where(
            Fdiff > 0,
            Fdiff * ((h[F1] - dem[F1]) / (h[F1] + 1)),
            Fdiff * ((h[F2] - dem[F2]) / (h[F2] + 1)),
        )
        Us[:, 1:-1, 1:-1] = (
            U[:, 1:-1, 1:-1] - tau * Ediff - tau * Fdiff - dt * S[:, 1:-1, 1:-1]
        )

    # ========================================== #
    def corrector(self):
        h = self.h
        dem = self.dem
        U = self.U
        Uss = self.Uss
        E = self.E
        F = self.F
        S = self.S
        alt = self.alt
        tau = self.tau
        dt = self.dt

        if alt == 0:
            E1 = np.s_[2:, 1:-1]
            E2 = np.s_[1:-1, 1:-1]
            F1 = np.s_[1:-1, 1:-1]
            F2 = np.s_[1:-1, 1:-1]
            Ediff = E[:, 2:, 1:-1] - E[:, 1:-1, 1:-1]
            Fdiff = F[:, 1:-1, 1:-1] - F[:, 1:-1, 1:-1]
        elif alt == 1:
            E1 = np.s_[1:-1, 1:-1]
            E2 = np.s_[1:-1, :-2]
            F1 = np.s_[1:-1, 1:-1]
            F2 = np.s_[1:-1, :-2]
            Ediff = E[:, 1:-1, 1:-1] - E[:, 1:-1, :-2]
            Fdiff = F[:, 1:-1, 1:-1] - F[:, 1:-1, :-2]
        elif alt == 2:
            E1 = np.s_[2:, 1:-1]
            E2 = np.s_[1:-1, 1:-1]
            F1 = np.s_[1:-1, 1:-1]
            F2 = np.s_[1:-1, :-2]
            Ediff = E[:, 2:, 1:-1] - E[:, 1:-1, 1:-1]
            Fdiff = F[:, 1:-1, 1:-1] - F[:, 1:-1, :-2]
        elif alt == 3:
            E1 = np.s_[1:-1, 1:-1]
            E2 = np.s_[:-2, 1:-1]
            F1 = np.s_[1:-1, 2:]
            F2 = np.s_[1:-1, 1:-1]
            Ediff = E[:, 1:-1, 1:-1] - E[:, :-2, 1:-1]
            Fdiff = F[:, 1:-1, 2:] - F[:, 1:-1, 1:-1]
        # -----------------------------------------------------------------------
        Ediff = np.where(
            Ediff > 0,
            Ediff * ((h[E1] - dem[E1]) / (h[E1] + 1)),
            Ediff * ((h[E2] - dem[E2]) / (h[E2] + 1)),
        )
        Fdiff = np.where(
            Fdiff > 0,
            Fdiff * ((h[F1] - dem[F1]) / (h[F1] + 1)),
            Fdiff * ((h[F2] - dem[F2]) / (h[F2] + 1)),
        )
        Uss[:, 1:-1, 1:-1] = (
            U[:, 1:-1, 1:-1] - tau * Ediff - tau * Fdiff - dt * S[:, 1:-1, 1:-1]
        )

    # ========================================== #
    def artificial_viscosity(self):

        """
        inputs: h,U,kappa,e,i,j,nx,ny
        """
        U = self.U
        h = self.h
        kappa = self.kappa
        e = self.e
        nx = self.nx
        ny = self.ny

        vis = np.zeros([nx, ny])
        epsilon_p = np.zeros([nx, ny])
        epsilon_m = np.zeros([nx, ny])
        epsilon_pp = np.zeros([nx, ny])
        epsilon_mm = np.zeros([nx, ny])

        vis_numerator = np.abs(h[2:, 1:-1] - 2 * h[1:-1, 1:-1] + h[:-2, 1:-1])

        vis_denominator = (
            e + np.abs(h[2:, 1:-1]) + 2 * np.abs(h[1:-1, 1:-1]) + np.abs(h[:-2, 1:-1])
        )
        vis[1:-1, 1:-1] = vis_numerator / vis_denominator

        epsilon_p[1:-1, 1:-1] = kappa * np.maximum(vis[2:, 1:-1], vis[1:-1, 1:-1])
        epsilon_m[1:-1, 1:-1] = kappa * np.maximum(vis[:-2, 1:-1], vis[1:-1, 1:-1])
        epsilon_pp[1:-1, 1:-1] = kappa * np.maximum(vis[1:-1, 2:], vis[1:-1, 1:-1])
        epsilon_mm[1:-1, 1:-1] = kappa * np.maximum(vis[1:-1, :-2], vis[1:-1, 1:-1])

        p = U[0, 2:, 1:-1] - U[0, 1:-1, 1:-1]
        m = U[0, 1:-1, 1:-1] - U[0, :-2, 1:-1]
        pp = U[0, 1:-1, 2:] - U[0, 1:-1, 1:-1]
        mm = U[0, 1:-1, 1:-1] - U[0, 1:-1, :-2]

        Utemp = U[0, 2:-2, 2:-2]

        Utemp = (
            Utemp
            + (
                epsilon_p[2:-2, 2:-2] * p[1:-1, 1:-1]
                - epsilon_m[2:-2, 2:-2] * (m[:-2, 1:-1])
                + epsilon_pp[2:-2, 2:-2] * (pp[1:-1, 1:-1])
                - epsilon_mm[2:-2, 2:-2] * (mm[1:-1, :-2])
            )
            / 2
        )

        U[1, 2:-2, 2:-2] = Utemp

    # ----------Solving for h and u,v and F,E and S-------
    def solve_for_h_u_v_F_E_and_s(self):

        g = self.g
        n = self.n
        h = self.h
        u = self.u
        v = self.v
        U = self.U
        F = self.F
        E = self.E
        S = self.S
        e = self.e

        h[1:-1, 1:-1] = U[0, 1:-1, 1:-1]
        u[1:-1, 1:-1] = U[1, 1:-1, 1:-1] / (U[0, 1:-1, 1:-1] + e)
        v[1:-1, 1:-1] = U[2, 1:-1, 1:-1] / (U[0, 1:-1, 1:-1] + e)

        E[0, :, :] = u * h
        E[1, :, :] = h * np.square(u) + g * np.square(h) / 2
        E[2, :, :] = u * v * h
        F[0, :, :] = v * h
        F[1, :, :] = u * v * h
        F[2, :, :] = h * np.square(v) + g * np.square(h) / 2

        S_numerator = g * h * np.abs(u) * u * np.square(n)
        S_denominator = e + np.power(np.maximum(h, np.zeros_like(h)), (1 / 3))
        S[1, :, :] = S_numerator / S_denominator

        S_numerator = g * h * np.abs(v) * self.v * np.square(n)
        S_denominator = e + np.power(np.maximum(h, np.zeros_like(h)), (1 / 3))
        S[2, :, :] = S_numerator / S_denominator

    # ========================================== #
    def set_synthetic_initial_head(self):
        self.dx = 1
        self.dy = 1
        self.xmax = 3
        self.ymax = 3
        self.x = np.arange(-self.xmax, self.xmax, self.dx)
        self.y = np.arange(-self.ymax, self.ymax, self.dy)
        self.nx = self.x.shape[0]
        self.ny = self.y.shape[0]
        self.h = np.zeros((self.nx, self.ny))
        for i in range(self.nx):
            for j in range(self.ny):
                self.h[i, j] = 10 * np.exp(
                    -(np.power(self.x[i], 2) + np.power(self.y[j], 2)) / 4
                )
        self.initial_wse = self.h
        self.dem = np.zeros_like(self.h) + np.random.rand(self.nx, self.ny)
        print(np.round(self.dem, 2))

    # ========================================== #
    def read_initial_head_from_file(self):
        # Will check if this is a data array to do the WSE
        da_hand_rt = xr.open_rasterio(self.data_dir + "hand_plus_rt.tif").astype(
            np.float64
        )
        da_dem = xr.open_rasterio(self.data_dir + "dem.tif").astype(np.float64)
        self.x = np.array(da_hand_rt.y.values)
        self.y = np.array(da_hand_rt.x.values)
        self.nx = self.x.shape[0]
        self.ny = self.y.shape[0]
        wse = np.where(da_hand_rt < 0.0, np.nan, da_hand_rt) + da_dem.fillna(0)

        self.initial_wse = wse.values[0, :, :]
        self.h = wse.values[0, :, :]

        # Set the initial depth, this will act as our boundary,
        self.initial_depth = da_hand_rt.values[0, :, :]
        self.initial_depth_boundary = np.where(self.initial_depth > 0, 1, 0)
        self.dem = da_dem.fillna(0).values[0, :, :]

        self.dx = 10
        self.dy = 10

    # ========================================== #
    def get_max_min_initial_head(self):
        self.max_h0 = np.max(self.h)
        self.min_h0 = np.min(self.h)

    # ========================================== #
    def set_U_F_E(self):

        self.U[0, :, :] = self.h
        self.F[2, :, :] = self.g * (np.square(self.h)) / 2
        self.E[1, :, :] = self.g * (np.square(self.h)) / 2
        self.htot = np.sum(self.h)
        self.alt = 0

    # -----------True U values---------%
    def get_true_U_values(self):

        U = self.U
        Us = self.Us
        Uss = self.Uss
        U[:, 1:-1, 1:-1] = (1 / 2) * (Us[:, 1:-1, 1:-1] + Uss[:, 1:-1, 1:-1])

    # ------------ Boundary Conditions---------------%
    def enforce_horizontal_boundary_conditions(self):
        self.U[0, :, :] = np.maximum(0, self.U[0, :, :])
        self.U[0, 0, :] = self.U[0, 1, :]
        self.U[0, :, 0] = self.U[0, :, 1]
        self.U[0, -1, :] = self.U[0, -2, :]
        self.U[0, :, -1] = self.U[0, :, -2]
        self.U[1, 0, :] = -self.U[1, 1, :]
        self.U[1, :, 0] = -self.U[1, :, 1]
        self.U[1, -1, :] = -self.U[1, -2, :]
        self.U[1, :, -1] = -self.U[1, :, -2]
        self.U[2, 0, :] = -self.U[2, 1, :]
        self.U[2, :, 0] = -self.U[1, :, 1]
        self.U[2, -1, :] = -self.U[2, -2, :]
        self.U[2, :, -1] = -self.U[1, :, -2]

    # ------------ Boundary Conditions---------------%
    def enforce_vertical_boundary_conditions(self):
        if self.vertical_boundary == "dem":
            self.h = np.maximum(self.h, self.dem)
        else:
            return

    # ========================================== #
    def set_alternative_gradient(self):
        self.alt = int(random.uniform(0, 4))

    # ========================================== #
    def update(self):
        self.t = self.t + self.dt
        self.count = self.count + 1
        self.nplot = self.nplot + 1

        self.set_alternative_gradient()

        # ---------------Predictor---------
        self.predictor()

        # ----------Corrector---------
        self.corrector()

        # -----------True U values---------%
        self.get_true_U_values()

        # -----------Artificial Viscosity---------------%
        self.artificial_viscosity()

        # ------------ Boundary Conditions---------------%
        self.enforce_vertical_boundary_conditions()
        self.enforce_horizontal_boundary_conditions()

        # ----------Solving for h and u,v and F,E and S-------
        self.solve_for_h_u_v_F_E_and_s()

    # ========================================== #
    def plot_head_at_time(self, save_plot_here=None):
        da = xr.DataArray(self.h).fillna(0)
        plt.figure()
        plt.contourf(
            da,
            levels=np.linspace(np.min(self.initial_wse), np.max(self.initial_wse), 20),
            # levels=np.linspace(self.min_h0, self.max_h0, 25),
            cmap="GnBu",
        )
        plt.colorbar()
        if save_plot_here is None:
            plt.show()
        else:
            plot_name = f"{save_plot_here}{np.round(self.t,3)}.tif"
            plt.savefig(plot_name)
            image = imageio.imread(plot_name)
            self.images.append(image)

        plt.close()

    # ========================================== #
    def _parse_config(self, cfg):
        for key, val in cfg.items():
            setattr(self, key, val)

    # ========================================== #
    def track_progress(self):
        if self.verbosity > 0:
            sv.plot_head_at_time(save_plot_here=sv.save_plots_here)
            print("time", sv.t)
            print("Mean h across domain", np.nanmean(sv.h))
            print("Sum h across domain", np.nansum(sv.h))


# ====================================================
# ====================================================
if __name__ == "__main__":

    sv = SAINT_VENANT()

    sv.initialize(config_file="sv_config.yml")

    sv.set_U_F_E()

    track_count = 0
    sv.track_progress()

    while sv.t < sv.tmax:

        if track_count == sv.ntplot:
            sv.track_progress()
            track_count = 0

        sv.update()

        if np.sum(sv.h) > 5 * sv.htot:
            sv.track_progress()
            break

        track_count += 1

    sv.track_progress()
    if sv.verbosity > 0:
        imageio.mimsave(sv.save_gif_here, sv.images)