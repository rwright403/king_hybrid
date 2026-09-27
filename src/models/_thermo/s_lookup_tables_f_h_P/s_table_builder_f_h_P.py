import numpy as np
import CoolProp.CoolProp as CP
import os


# ============================================================
# OUTPUT FILE
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

OUT_FILE = os.path.join(
    SCRIPT_DIR,
    "n2o_s_table_f_h_P.npz"
)


# ============================================================
# BUILD TABLE
# ============================================================

def build_n2o_s_table(
    T_min=183.0,          # K
    T_max=309.0,          # K
    P_min=0.7e6,          # Pa
    P_max=7.0e6,          # Pa
    nT_envelope=200,
    nh=500,
    nP=200,
    out_file=OUT_FILE,
):
    """
    Build an N2O entropy lookup table:

        s = f(h, P)

    using CoolProp.

    The T-P envelope is used to determine a COMMON enthalpy
    range that exists across the full pressure range.

    A rectangular h-P lookup table is then generated.

    Units:
        T : K
        P : Pa
        h : J/kg
        s : J/(kg K)

    Saved variables:
        h_grid : enthalpy grid [J/kg]
        P_grid : pressure grid [Pa]
        s_tab  : entropy table [J/(kg K)]

    Table indexing:
        s_tab[i, j] = s(h_grid[i], P_grid[j])
    """

    fluid = "N2O"


    # ========================================================
    # PRESSURE GRID
    # ========================================================

    P_grid = np.linspace(
        P_min,
        P_max,
        nP,
        dtype=np.float64
    )


    # ========================================================
    # TEMPERATURE ENVELOPE
    # ========================================================

    T_envelope = np.linspace(
        T_min,
        T_max,
        nT_envelope,
        dtype=np.float64
    )


    # ========================================================
    # FIND ENTHALPY RANGE AT EACH PRESSURE
    # ========================================================

    print("Determining enthalpy envelope...")

    h_min_by_P = []
    h_max_by_P = []

    envelope_failures = 0


    for j, P in enumerate(P_grid):

        h_at_P = []

        for T in T_envelope:

            try:

                h = CP.PropsSI(
                    "H",
                    "T", T,
                    "P", P,
                    fluid
                )

                if np.isfinite(h):
                    h_at_P.append(h)

            except Exception:
                envelope_failures += 1


        if len(h_at_P) == 0:

            raise RuntimeError(
                f"No valid enthalpy states found at "
                f"P = {P / 1e6:.4f} MPa"
            )


        h_min_by_P.append(np.min(h_at_P))
        h_max_by_P.append(np.max(h_at_P))


        print(
            f"\rEnvelope progress: {j + 1}/{nP} "
            f"({100 * (j + 1) / nP:.1f}%)",
            end=""
        )


    print("\n")


    h_min_by_P = np.asarray(
        h_min_by_P,
        dtype=np.float64
    )

    h_max_by_P = np.asarray(
        h_max_by_P,
        dtype=np.float64
    )


    # ========================================================
    # FIND COMMON ENTHALPY RANGE
    # ========================================================

    # For a rectangular h-P table, we want an enthalpy range
    # that exists at EVERY pressure.
    #
    # Therefore:
    #
    #   h_min = highest lower bound
    #   h_max = lowest upper bound

    h_min = np.max(h_min_by_P)
    h_max = np.min(h_max_by_P)


    if h_min >= h_max:

        raise RuntimeError(
            "No common enthalpy range exists across the "
            "requested pressure range.\n"
            f"h_min = {h_min:.3f} J/kg\n"
            f"h_max = {h_max:.3f} J/kg"
        )


    print("Thermodynamic envelope:")

    print(
        f"    Temperature : "
        f"{T_min:.2f} - {T_max:.2f} K"
    )

    print(
        f"    Pressure    : "
        f"{P_min / 1e6:.3f} - "
        f"{P_max / 1e6:.3f} MPa"
    )

    print(
        f"    Enthalpy    : "
        f"{h_min:.3f} - "
        f"{h_max:.3f} J/kg"
    )

    print(
        f"                  "
        f"({h_min / 1e3:.3f} - "
        f"{h_max / 1e3:.3f} kJ/kg)"
    )

    print(
        f"    Envelope CoolProp failures: "
        f"{envelope_failures}"
    )

    print()


    # ========================================================
    # ENTHALPY GRID
    # ========================================================

    h_grid = np.linspace(
        h_min,
        h_max,
        nh,
        dtype=np.float64
    )


    # ========================================================
    # BUILD ENTROPY TABLE
    # ========================================================

    # Table orientation:
    #
    #                   Pressure
    #                     j ->
    #
    # Enthalpy     [ s  s  s  s ]
    #    i         [ s  s  s  s ]
    #    |         [ s  s  s  s ]
    #    v         [ s  s  s  s ]
    #
    #
    # s_tab[i, j]
    #     = s(h_grid[i], P_grid[j])

    s_tab = np.full(
        (nh, nP),
        np.nan,
        dtype=np.float64
    )


    print("Building entropy table...")


    table_failures = 0

    # Only print the first few failures
    max_failure_prints = 10


    for i, h in enumerate(h_grid):

        for j, P in enumerate(P_grid):

            try:

                s = CP.PropsSI(
                    "S",
                    "H", h,
                    "P", P,
                    fluid
                )

                if np.isfinite(s):

                    s_tab[i, j] = s

                else:

                    table_failures += 1

                    if table_failures <= max_failure_prints:

                        print(
                            f"\nNon-finite result:"
                            f" h={h / 1e3:.3f} kJ/kg,"
                            f" P={P / 1e6:.4f} MPa"
                        )


            except Exception as e:

                table_failures += 1

                s_tab[i, j] = np.nan

                if table_failures <= max_failure_prints:

                    print()

                    print(
                        f"FAILED:"
                        f" h={h / 1e3:.3f} kJ/kg,"
                        f" P={P / 1e6:.4f} MPa"
                    )

                    print(
                        f"CoolProp error: {e}"
                    )


        print(
            f"\rProgress: {i + 1}/{nh} "
            f"({100 * (i + 1) / nh:.1f}%)",
            end=""
        )


    print("\n")


    # ========================================================
    # TABLE VALIDATION
    # ========================================================

    total_points = s_tab.size

    valid_mask = np.isfinite(s_tab)

    valid_points = np.count_nonzero(
        valid_mask
    )

    invalid_points = (
        total_points - valid_points
    )

    invalid_fraction = (
        invalid_points
        / total_points
        * 100.0
    )


    print("Table complete:")

    print(
        f"    Shape          : {s_tab.shape}"
    )

    print(
        f"    Total points   : {total_points}"
    )

    print(
        f"    Valid points   : {valid_points}"
    )

    print(
        f"    Invalid points : {invalid_points}"
    )

    print(
        f"    Invalid        : {invalid_fraction:.3f} %"
    )

    print()


    # ========================================================
    # REPORT INVALID LOCATIONS
    # ========================================================

    if invalid_points > 0:

        bad = np.argwhere(
            ~valid_mask
        )

        print("First invalid table states:")

        for i, j in bad[:10]:

            print(
                f"    h = "
                f"{h_grid[i] / 1e3:.3f} kJ/kg, "
                f"P = "
                f"{P_grid[j] / 1e6:.4f} MPa"
            )

        print()


    # ========================================================
    # VALIDITY BY PRESSURE
    # ========================================================

    print("Validity by pressure:")

    for j, P in enumerate(P_grid):

        valid_at_P = np.count_nonzero(
            np.isfinite(s_tab[:, j])
        )

        # Only print pressure columns containing invalid states
        if valid_at_P != nh:

            print(
                f"    P = {P / 1e6:.4f} MPa: "
                f"{valid_at_P}/{nh} valid"
            )


    if invalid_points == 0:

        print(
            "    All pressure columns are fully valid."
        )

    print()


    # ========================================================
    # SAVE TABLE
    # ========================================================

    np.savez(
        out_file,
        h_grid=h_grid,
        P_grid=P_grid,
        s_tab=s_tab,
    )


    print("Saved to:")
    print(f"    {out_file}")


    # ========================================================
    # RETURN
    # ========================================================

    return h_grid, P_grid, s_tab


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    h_grid, P_grid, s_tab = build_n2o_s_table()