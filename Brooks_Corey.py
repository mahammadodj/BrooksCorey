import numpy as np
import matplotlib.pyplot as plt


class ReservoirModel:
    def __init__(self, ro_aqu, ro_oleic, WOC, GOC,
                 M_gas=19, Z=1.2, R=10.73, T=175,
                 sigma_ow=20, sigma_og=30, sigma_gw=50,
                 g=1 / 144, Swirr=0.2, _lambda=2, Pe_ow=3.5):
        """
        ro_aqu, ro_oleic : water / oil density [lb/ft3]
        WOC, GOC         : contact depths [ft] (TVD from the pressure datum)
        M_gas            : gas molecular weight [lb/lbmol]
        Z                : gas compressibility factor
        R                : gas constant [psia*ft3/(lbmol*R)]
        T                : reservoir temperature [F]
        sigma_*          : interfacial tensions [dyne/cm]
        g                : unit conversion, lb/ft3 -> psi/ft
        Swirr            : irreducible water saturation
        _lambda          : Brooks-Corey pore size distribution index
        Pe_ow            : oil-water entry pressure [psi]
        """
        self.ro_aqu = ro_aqu
        self.ro_oleic = ro_oleic
        self.M_gas = M_gas
        self.Z = Z
        self.R = R
        self.T = T
        self.sigma_ow = sigma_ow
        self.sigma_og = sigma_og
        self.sigma_gw = sigma_gw
        self.g = g
        self.Swirr = Swirr
        self._lambda = _lambda
        self.Pe_ow = Pe_ow
        self.WOC = WOC
        self.GOC = GOC

    # ------------------------------------------------------------------
    # Entry pressures
    # ------------------------------------------------------------------
    def calculate_entry_pressures(self):
        Pe_ow = self.Pe_ow
        Pe_og = Pe_ow * self.sigma_og / self.sigma_ow
        Pe_wg = Pe_ow * self.sigma_gw / self.sigma_ow
        return Pe_ow, Pe_og, Pe_wg

    # ------------------------------------------------------------------
    # Pressures at the contacts
    # ------------------------------------------------------------------
    def calculate_p_contacts(self):
        Pe_ow, Pe_og, _ = self.calculate_entry_pressures()

        # Water pressure at WOC, oil pressure is shifted by the entry pressure
        Pw_WOC = 14.7 + self.ro_aqu * self.g * self.WOC
        Po_WOC = Pw_WOC + Pe_ow

        # Oil pressure at GOC (water is not needed here)
        Po_GOC = Po_WOC + self.ro_oleic * self.g * (self.GOC - self.WOC)

        # Gas density evaluated at the actual oil pressure at the GOC
        # (the original rounding up to the next 100 psi was removed)
        T_rankine = self.T + 459.67
        ro_g = (Po_GOC * self.M_gas) / (self.Z * self.R * T_rankine)

        # Gas pressure at GOC
        Pg_GOC = Po_GOC + Pe_og

        return Pw_WOC, Po_WOC, Po_GOC, Pg_GOC, ro_g

    def calculate_FWL(self):
        """
        FWL is the depth where oil and water pressures are equal (Pc_ow = 0).
        At the WOC, Po - Pw = Pe_ow, and the difference shrinks with depth
        at (ro_aqu - ro_oleic) * g psi/ft, so:
            FWL = WOC + Pe_ow / ((ro_aqu - ro_oleic) * g)
        """
        grad_diff = (self.ro_aqu - self.ro_oleic) * self.g
        if grad_diff <= 0:
            raise ValueError("Water must be denser than oil to define an FWL.")
        return self.WOC + self.Pe_ow / grad_diff

    def calculate_FGL(self):
        """
        Free gas level: depth where gas and oil pressures are equal (Pc_go = 0).
        At the GOC, Pg - Po = Pe_og, and the difference shrinks with depth
        at (ro_oleic - ro_g) * g psi/ft.
        """
        _, _, _, _, ro_g = self.calculate_p_contacts()
        _, Pe_og, _ = self.calculate_entry_pressures()
        grad_diff = (self.ro_oleic - ro_g) * self.g
        if grad_diff <= 0:
            raise ValueError("Oil must be denser than gas to define an FGL.")
        return self.GOC + Pe_og / grad_diff

    # ------------------------------------------------------------------
    # Phase pressures at any depth
    # ------------------------------------------------------------------
    def calculate_pressures(self, Depth):
        Pw_WOC, Po_WOC, _, Pg_GOC, ro_g = self.calculate_p_contacts()
        Po = Po_WOC + self.ro_oleic * self.g * (Depth - self.WOC)
        Pw = Pw_WOC + self.ro_aqu * self.g * (Depth - self.WOC)
        Pg = Pg_GOC + ro_g * self.g * (Depth - self.GOC)
        return Po, Pw, Pg

    # ------------------------------------------------------------------
    # Brooks-Corey effective saturation
    # ------------------------------------------------------------------
    def _brooks_corey(self, Pc, Pe):
        """Wetting-phase saturation for capillary pressure Pc and entry pressure Pe."""
        if Pc <= Pe:  # below entry pressure -> fully wetting phase
            return 1.0
        return self.Swirr + (1 - self.Swirr) * (Pc / Pe) ** (-self._lambda)

   
    # Saturations at any depth
    def calculate_saturations(self, Depth, Po, Pw, Pg):
        Pe_ow, Pe_og, _ = self.calculate_entry_pressures()

        Sw = self._brooks_corey(Po - Pw, Pe_ow)   # water saturation (oil-water Pc)
        Sl = self._brooks_corey(Pg - Po, Pe_og)   # total liquid saturation (gas-oil Pc)
        Sl = max(Sl, Sw)                          # liquid can never be less than water

        So = Sl - Sw
        Sg = 1 - Sl
        return Sw, So, Sg



if __name__ == "__main__":
    WOC = 7475
    GOC = 7450

    reservoir = ReservoirModel(
        WOC=WOC,
        GOC=GOC,
        ro_aqu=66.96,
        ro_oleic=53
    )

    FWL = reservoir.calculate_FWL()
    FGL = reservoir.calculate_FGL()

    print(f"Free Water Level (FWL): {FWL:.1f} ft")
    print(f"Free Gas Level  (FGL): {FGL:.1f} ft")
    print(f"FWL - WOC = {FWL - WOC:.1f} ft")

    depths = np.arange(7300, 7550)

    # Saturation lists
    Sw_, So_, Sg_ = [], [], []

    # Pressure lists and corresponding depths
    Pw_, Po_, Pg_ = [], [], []
    depths_water = []
    depths_oil = []
    depths_gas = []

    for depth in depths:
        Po, Pw, Pg = reservoir.calculate_pressures(depth)

        Sw, So, Sg = reservoir.calculate_saturations(
            Depth=depth, Po=Po, Pw=Pw, Pg=Pg
        )

        # Store saturations
        Sw_.append(Sw)
        So_.append(So)
        Sg_.append(Sg)

        # Store water pressure and its corresponding depth
        Pw_.append(Pw)
        depths_water.append(depth)

        # Store oil pressure only up to FWL
        if depth <= FWL:
            Po_.append(Po)
            depths_oil.append(depth)

        # Store gas pressure only up to FGL
        if depth <= FGL:
            Pg_.append(Pg)
            depths_gas.append(depth)

    # Check saturation balance
    assert np.allclose(
        np.array(Sw_) + np.array(So_) + np.array(Sg_),
        1.0
    ), "Saturations do not sum to 1"
    

    # -------------------------------
    # Plot saturations
    # -------------------------------
    plt.figure(figsize=(6, 8))

    plt.plot(Sw_, depths, color="blue", label="Water Saturation")
    plt.plot(So_, depths, color="red", label="Oil Saturation")
    plt.plot(Sg_, depths, color="green", label="Gas Saturation")

    plt.ylim(7550, 7300)
    plt.xlim(0, 1)

    for level, label in [
        (GOC, "GOC"),
        (WOC, "WOC"),
        (FGL, "FGL"),
        (FWL, "FWL")
    ]:
        plt.axhline(y=level, linestyle="--", color="black", alpha=0.5)
        plt.text(0.9, level - 4, label, fontsize=10)

    plt.xlabel("Phase Saturation")
    plt.ylabel("Depth [ft]")
    plt.title("Fluid Saturation vs Depth")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
    

    # -------------------------------
    # Plot pressures
    # -------------------------------
    plt.figure(figsize=(6, 8))

    plt.plot(Pw_, depths_water, color="blue", label="Water Pressure")
    plt.plot(Po_, depths_oil, color="red", label="Oil Pressure")
    plt.plot(Pg_, depths_gas, color="green", label="Gas Pressure")

    plt.ylim(7550, 7300)

    for level, label in [
        (GOC, "GOC"),
        (WOC, "WOC"),
        (FGL, "FGL"),
        (FWL, "FWL")
    ]:
        plt.axhline(y=level, linestyle="--", color="black", alpha=0.4)
        plt.text(
            plt.xlim()[1], level - 4, label,
            fontsize=9, ha="right"
        )

    plt.xlabel("Phase Pressure [psi]")
    plt.ylabel("Depth [ft]")
    plt.title("Fluid Pressure vs Depth")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
