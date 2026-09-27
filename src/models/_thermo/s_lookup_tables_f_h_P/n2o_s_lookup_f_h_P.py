import numpy as np
import CoolProp.CoolProp as CP


class N2O_s_f_h_P_Table:
    def __init__(
        self,
        path="src/models/_thermo/s_lookup_tables_f_h_P/n2o_s_table_f_h_P.npz"
    ):
        d = np.load(path)

        self.h_grid = d["h_grid"]
        self.P_grid = d["P_grid"]
        self.s_tab = d["s_tab"]


    def lookup(self, h, P):
        h_grid = self.h_grid
        P_grid = self.P_grid
        F = self.s_tab

        i = np.searchsorted(h_grid, h) - 1
        j = np.searchsorted(P_grid, P) - 1

        i = max(0, min(i, len(h_grid) - 2))
        j = max(0, min(j, len(P_grid) - 2))

        H1, H2 = h_grid[i], h_grid[i + 1]
        P1, P2 = P_grid[j], P_grid[j + 1]

        Q11 = F[i,     j]
        Q21 = F[i + 1, j]
        Q12 = F[i,     j + 1]
        Q22 = F[i + 1, j + 1]

        # Fallback if any NaN
        if (
            np.isnan(Q11)
            or np.isnan(Q21)
            or np.isnan(Q12)
            or np.isnan(Q22)
        ):
            print(
                f"Caution: entropy lookup returned NaN at "
                f"h={h:.0f} J/kg, P={P:.0f} Pa. "
                f"Using CoolProp directly."
            )

            return CP.PropsSI(
                "S",
                "H", h,
                "P", P,
                "N2O"
            )

        # Bilinear interpolation
        a = (h - H1) / (H2 - H1)
        b = (P - P1) / (P2 - P1)

        return (
            (1 - a) * (1 - b) * Q11
            + a * (1 - b) * Q21
            + (1 - a) * b * Q12
            + a * b * Q22
        )


import matplotlib.pyplot as plt


def plot_lookup_error(table, nh=200, nP=200):
    """
    Plot percent error between the lookup-table interpolation
    and direct CoolProp values over the full h-P table domain.

    x-axis : enthalpy [kJ/kg]
    y-axis : pressure [MPa]
    color  : absolute percent error [%]
    """

    # Test points across the table domain.
    # Offset slightly from the exact edges.
    h_test = np.linspace(
        table.h_grid[0],
        table.h_grid[-1],
        nh
    )

    P_test = np.linspace(
        table.P_grid[0],
        table.P_grid[-1],
        nP
    )

    error = np.full(
        (nP, nh),
        np.nan,
        dtype=np.float64
    )

    print("Calculating interpolation error...")

    for j, P in enumerate(P_test):

        for i, h in enumerate(h_test):

            try:
                # Lookup table
                s_interp = table.lookup(h, P)

                # Direct CoolProp
                s_CP = CP.PropsSI(
                    "S",
                    "H", h,
                    "P", P,
                    "N2O"
                )

                # Absolute percent error
                error[j, i] = (
                    abs(s_interp - s_CP)
                    / abs(s_CP)
                    * 100.0
                )

            except Exception:
                error[j, i] = np.nan

        print(
            f"\rProgress: {j + 1}/{nP} "
            f"({100 * (j + 1) / nP:.1f}%)",
            end=""
        )

    print("\n")


    # ========================================================
    # ERROR STATISTICS
    # ========================================================

    valid_error = error[np.isfinite(error)]

    print("Interpolation error:")
    print(f"    Mean : {np.mean(valid_error):.6f} %")
    print(f"    95th : {np.percentile(valid_error, 95):.6f} %")
    print(f"    Max  : {np.max(valid_error):.6f} %")


    # ========================================================
    # PLOT
    # ========================================================

    plt.figure(figsize=(10, 7))

    mesh = plt.pcolormesh(
        h_test / 1e3,      # kJ/kg
        P_test / 1e6,      # MPa
        error,
        shading="auto"
    )

    cbar = plt.colorbar(mesh)
    cbar.set_label("Absolute Error [%]")

    plt.xlabel("Enthalpy [kJ/kg]")
    plt.ylabel("Pressure [MPa]")

    plt.title(
        "N₂O Entropy Lookup Table Error\n"
        "Bilinear Interpolation vs CoolProp"
    )

    plt.tight_layout()
    plt.show()


# ============================================================
# RUN VALIDATION
# ============================================================

if __name__ == "__main__":

    table = N2O_S_f_h_P_GasTable()

    plot_lookup_error(
        table,
        nh=200,
        nP=200
    )