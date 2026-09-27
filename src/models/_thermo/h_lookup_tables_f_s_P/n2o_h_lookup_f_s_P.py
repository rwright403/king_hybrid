import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt


# ============================================================
# LOOKUP TABLE
# ============================================================

class N2O_h_f_s_P_Table:
    def __init__(
        self,
        path="src/models/_thermo/h_lookup_tables_f_s_P/n2o_h_table_f_s_P.npz"
    ):
        d = np.load(path)

        self.s_grid = d["s_grid"]
        self.P_grid = d["P_grid"]
        self.h_tab = d["h_tab"]


    def lookup(self, s, P):
        s_grid = self.s_grid
        P_grid = self.P_grid
        F = self.h_tab

        i = np.searchsorted(s_grid, s) - 1
        j = np.searchsorted(P_grid, P) - 1

        i = max(0, min(i, len(s_grid) - 2))
        j = max(0, min(j, len(P_grid) - 2))

        S1, S2 = s_grid[i], s_grid[i + 1]
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
                f"Caution: enthalpy lookup returned NaN at "
                f"s={s:.2f} J/(kg K), P={P:.0f} Pa. "
                f"Using CoolProp directly."
            )

            return CP.PropsSI(
                "H",
                "S", s,
                "P", P,
                "N2O"
            )

        # Bilinear interpolation
        a = (s - S1) / (S2 - S1)
        b = (P - P1) / (P2 - P1)

        return (
            (1 - a) * (1 - b) * Q11
            + a * (1 - b) * Q21
            + (1 - a) * b * Q12
            + a * b * Q22
        )


# ============================================================
# ERROR PLOT
# ============================================================

def plot_lookup_error(table, ns=200, nP=200):
    """
    Plot percent error between the lookup-table interpolation
    and direct CoolProp values over the full s-P table domain.

    x-axis : entropy [kJ/(kg K)]
    y-axis : pressure [MPa]
    color  : absolute percent error [%]
    """

    # Test points across the table domain
    s_test = np.linspace(
        table.s_grid[0],
        table.s_grid[-1],
        ns
    )

    P_test = np.linspace(
        table.P_grid[0],
        table.P_grid[-1],
        nP
    )

    error = np.full(
        (nP, ns),
        np.nan,
        dtype=np.float64
    )

    print("Calculating interpolation error...")

    for j, P in enumerate(P_test):

        for i, s in enumerate(s_test):

            try:
                # Lookup table
                h_interp = table.lookup(s, P)

                # Direct CoolProp
                h_CP = CP.PropsSI(
                    "H",
                    "S", s,
                    "P", P,
                    "N2O"
                )

                # Absolute percent error
                error[j, i] = (
                    abs(h_interp - h_CP)
                    / abs(h_CP)
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
        s_test / 1e3,      # kJ/(kg K)
        P_test / 1e6,      # MPa
        error,
        shading="auto"
    )

    cbar = plt.colorbar(mesh)
    cbar.set_label("Absolute Error [%]")

    plt.xlabel("Entropy [kJ/(kg K)]")
    plt.ylabel("Pressure [MPa]")

    plt.title(
        "N₂O Enthalpy Lookup Table Error\n"
        "Bilinear Interpolation vs CoolProp"
    )

    plt.tight_layout()
    plt.show()


# ============================================================
# RUN VALIDATION
# ============================================================

if __name__ == "__main__":

    table = N2O_h_f_s_P_Table()

    plot_lookup_error(
        table,
        ns=200,
        nP=200
    )