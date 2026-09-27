import numpy as np
import CoolProp.CoolProp as CP
import os


# ============================================================
# OUTPUT FILE
# ============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

OUT_FILE = os.path.join(
    SCRIPT_DIR,
    "n2o_rho_table_f_s_P.npz"
)


# ============================================================
# BUILD TABLE
# ============================================================

def build_n2o_rho_table(
    T_min=183.0,          # K
    T_max=309.0,          # K
    P_min=0.7e6,          # Pa
    P_max=7.0e6,          # Pa
    nT_envelope=200,
    ns=500,
    nP=200,
    out_file=OUT_FILE,
):
    """
    Build an N2O density lookup table:

        rho = f(s, P)

    using CoolProp.

    The T-P envelope is used to determine a COMMON entropy
    range that exists across the full pressure range.

    A rectangular s-P lookup table is then generated.

    Units:
        T   : K
        P   : Pa
        s   : J/(kg K)
        rho : kg/m^3

    Saved variables:
        s_grid  : entropy grid [J/(kg K)]
        P_grid  : pressure grid [Pa]
        rho_tab : density table [kg/m^3]

    Table indexing:
        rho_tab[i, j] = rho(s_grid[i], P_grid[j])
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
    # FIND ENTROPY RANGE AT EACH PRESSURE
    # ========================================================

    print("Determining entropy envelope...")

    s_min_by_P = []
    s_max_by_P = []

    envelope_failures = 0


    for j, P in enumerate(P_grid):

        s_at_P = []

        for T in T_envelope:

            try:

                s = CP.PropsSI(
                    "S",
                    "T", T,
                    "P", P,
                    fluid
                )

                if np.isfinite(s):
                    s_at_P.append(s)

            except Exception:
                envelope_failures += 1


        if len(s_at_P) == 0:

            raise RuntimeError(
                f"No valid entropy states found at "
                f"P = {P / 1e6:.4f} MPa"
            )


        s_min_by_P.append(np.min(s_at_P))
        s_max_by_P.append(np.max(s_at_P))


        print(
            f"\rEnvelope progress: {j + 1}/{nP} "
            f"({100 * (j + 1) / nP:.1f}%)",
            end=""
        )


    print("\n")


    s_min_by_P = np.asarray(
        s_min_by_P,
        dtype=np.float64
    )

    s_max_by_P = np.asarray(
        s_max_by_P,
        dtype=np.float64
    )


    # ========================================================
    # FIND COMMON ENTROPY RANGE
    # ========================================================

    # We want a rectangular s-P domain that is valid across
    # the entire pressure range.
    #
    # Therefore:
    #
    #   s_min = highest lower bound
    #   s_max = lowest upper bound

    s_min = np.max(s_min_by_P)
    s_max = np.min(s_max_by_P)


    if s_min >= s_max:

        raise RuntimeError(
            "No common entropy range exists across the "
            "requested pressure range.\n"
            f"s_min = {s_min:.3f} J/(kg K)\n"
            f"s_max = {s_max:.3f} J/(kg K)"
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
        f"    Entropy     : "
        f"{s_min:.3f} - "
        f"{s_max:.3f} J/(kg K)"
    )

    print(
        f"                  "
        f"({s_min / 1e3:.6f} - "
        f"{s_max / 1e3:.6f} kJ/(kg K))"
    )

    print(
        f"    Envelope CoolProp failures: "
        f"{envelope_failures}"
    )

    print()


    # ========================================================
    # ENTROPY GRID
    # ========================================================

    s_grid = np.linspace(
        s_min,
        s_max,
        ns,
        dtype=np.float64
    )


    # ========================================================
    # BUILD DENSITY TABLE
    # ========================================================

    # Table orientation:
    #
    #                   Pressure
    #                     j ->
    #
    # Entropy      [ rho rho rho rho ]
    #    i         [ rho rho rho rho ]
    #    |         [ rho rho rho rho ]
    #    v         [ rho rho rho rho ]
    #
    #
    # rho_tab[i, j]
    #     = rho(s_grid[i], P_grid[j])

    rho_tab = np.full(
        (ns, nP),
        np.nan,
        dtype=np.float64
    )


    print("Building density table...")


    table_failures = 0

    # Only print the first few failures
    max_failure_prints = 10


    for i, s in enumerate(s_grid):

        for j, P in enumerate(P_grid):

            try:

                rho = CP.PropsSI(
                    "D",
                    "S", s,
                    "P", P,
                    fluid
                )

                if np.isfinite(rho):

                    rho_tab[i, j] = rho

                else:

                    table_failures += 1

                    if table_failures <= max_failure_prints:

                        print(
                            f"\nNon-finite result:"
                            f" s={s:.3f} J/(kg K),"
                            f" P={P / 1e6:.4f} MPa"
                        )


            except Exception as e:

                table_failures += 1

                rho_tab[i, j] = np.nan

                if table_failures <= max_failure_prints:

                    print()

                    print(
                        f"FAILED:"
                        f" s={s:.3f} J/(kg K),"
                        f" P={P / 1e6:.4f} MPa"
                    )

                    print(
                        f"CoolProp error: {e}"
                    )


        print(
            f"\rProgress: {i + 1}/{ns} "
            f"({100 * (i + 1) / ns:.1f}%)",
            end=""
        )


    print("\n")


    # ========================================================
    # TABLE VALIDATION
    # ========================================================

    total_points = rho_tab.size

    valid_mask = np.isfinite(rho_tab)

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
        f"    Shape          : {rho_tab.shape}"
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
                f"    s = "
                f"{s_grid[i]:.3f} J/(kg K), "
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
            np.isfinite(rho_tab[:, j])
        )

        if valid_at_P != ns:

            print(
                f"    P = {P / 1e6:.4f} MPa: "
                f"{valid_at_P}/{ns} valid"
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
        s_grid=s_grid,
        P_grid=P_grid,
        rho_tab=rho_tab,
    )


    print("Saved to:")
    print(f"    {out_file}")


    # ========================================================
    # RETURN
    # ========================================================

    return s_grid, P_grid, rho_tab


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    s_grid, P_grid, rho_tab = build_n2o_rho_table()