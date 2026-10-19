import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


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

    def calculate_entry_pressures(self):
        Pe_ow = self.Pe_ow
        Pe_og = Pe_ow * self.sigma_og / self.sigma_ow
        Pe_wg = Pe_ow * self.sigma_gw / self.sigma_ow
        return Pe_ow, Pe_og, Pe_wg

    
    def calculate_p_contacts(self):
        Pe_ow, Pe_og, _ = self.calculate_entry_pressures()

        # Water pressure at WOC, oil pressure is shifted by the entry pressure
        Pw_WOC = 14.7 + self.ro_aqu * self.g * self.WOC
        Po_WOC = Pw_WOC + Pe_ow

        # Oil pressure at GOC (water is not needed here)
        Po_GOC = Po_WOC + self.ro_oleic * self.g * (self.GOC - self.WOC)

        # Gas density evaluated at the actual oil pressure at the GOC
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


    def calculate_saturations(self, Depth, Po, Pw, Pg):
        Pe_ow, Pe_og, _ = self.calculate_entry_pressures()

        Sw = self._brooks_corey(Po - Pw, Pe_ow)   
        Sl = self._brooks_corey(Pg - Po, Pe_og)   
        Sl = max(Sl, Sw)

        So = Sl - Sw
        Sg = 1 - Sl
        return Sw, So, Sg



if __name__ == "__main__":
    WOC = 7475
    GOC = 7450

    reservoir = ReservoirModel(WOC=WOC, GOC=GOC, ro_aqu=66.96, ro_oleic=53)

    FWL = reservoir.calculate_FWL()
    FGL = reservoir.calculate_FGL()
    print(f"Free Water Level (FWL): {FWL:.1f} ft")
    print(f"Free Gas Level  (FGL): {FGL:.1f} ft")
    print(f"FWL - WOC = {FWL - WOC:.1f} ft")

    depths = np.arange(7300, 7550)

    Sw_, So_, Sg_ = [], [], []
    Pw_, Po_, Pg_ = [], [], []
    for depth in depths:
        Po, Pw, Pg = reservoir.calculate_pressures(depth)
        Sw, So, Sg = reservoir.calculate_saturations(Depth=depth, Po=Po, Pw=Pw, Pg=Pg)
        Sw_.append(Sw); So_.append(So); Sg_.append(Sg)
        Pw_.append(Pw); Po_.append(Po); Pg_.append(Pg)

    Sw_, So_, Sg_ = map(np.array, (Sw_, So_, Sg_))
    Pw_, Po_, Pg_ = map(np.array, (Pw_, Po_, Pg_))

    assert np.allclose(Sw_ + So_ + Sg_, 1.0), "Saturations do not sum to 1"

    n_oil = int(np.sum(depths <= FWL))
    n_gas = int(np.sum(depths <= FGL))

    fig, (ax_s, ax_p) = plt.subplots(1, 2, figsize=(12, 8), sharey=True)

    ax_s.set_xlim(0, 1)
    ax_s.set_ylim(depths.max(), depths.min())
    ax_s.set_xlabel("Phase Saturation")
    ax_s.set_ylabel("Depth [ft]")
    ax_s.set_title("Fluid Saturation vs Depth")

    pmin = min(Pw_.min(), Po_[:n_oil].min(), Pg_[:n_gas].min()) - 5
    pmax = max(Pw_.max(), Po_[:n_oil].max(), Pg_[:n_gas].max()) + 5
    ax_p.set_xlim(pmin, pmax)
    ax_p.set_xlabel("Phase Pressure [psi]")
    ax_p.set_title("Fluid Pressure vs Depth")

    # Contact / free-level lines
    levels = [(GOC, "GOC"), (WOC, "WOC"), (FGL, "FGL"), (FWL, "FWL")]
    for a, x_text, ha in [(ax_s, 0.98, "right"), (ax_p, pmax - 1, "right")]:
        for level, name in levels:
            a.axhline(level, linestyle="--", color="black", alpha=0.4)
            a.text(x_text, level - 3, name, fontsize=9, ha=ha)
        a.grid(alpha=0.3)

    # Animated curves (start empty)
    (l_sw,) = ax_s.plot([], [], color="blue", label="Water Saturation")
    (l_so,) = ax_s.plot([], [], color="red", label="Oil Saturation")
    (l_sg,) = ax_s.plot([], [], color="green", label="Gas Saturation")
    (l_pw,) = ax_p.plot([], [], color="blue", label="Water Pressure")
    (l_po,) = ax_p.plot([], [], color="red", label="Oil Pressure")
    (l_pg,) = ax_p.plot([], [], color="green", label="Gas Pressure")

    # Moving markers
    (d_sw,) = ax_s.plot([], [], "o", color="blue")
    (d_so,) = ax_s.plot([], [], "o", color="red")
    (d_sg,) = ax_s.plot([], [], "o", color="green")
    (d_pw,) = ax_p.plot([], [], "o", color="blue")
    (d_po,) = ax_p.plot([], [], "o", color="red")
    (d_pg,) = ax_p.plot([], [], "o", color="green")

    # Moving value labels next to the saturation dots (stacked so they never overlap)
    def make_label(color, y_off):
        return ax_s.annotate("", xy=(0, 0), xytext=(10, y_off), textcoords="offset points",
                             color=color, fontsize=10, fontweight="bold", va="center",
                             bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                       edgecolor=color, alpha=0.9))
    lab_sw = make_label("blue", -18)
    lab_so = make_label("red", 0)
    lab_sg = make_label("green", 18)
    sat_labels = [(lab_sw, Sw_, "Sw"), (lab_so, So_, "So"), (lab_sg, Sg_, "Sg")]

    # Depth cursors
    cur_s = ax_s.axhline(depths[0], color="grey", linewidth=1.5)
    cur_p = ax_p.axhline(depths[0], color="grey", linewidth=1.5)

    # Read-out boxes
    box = dict(boxstyle="round", facecolor="white", alpha=0.85)
    txt_s = ax_s.text(0.03, 0.03, "", transform=ax_s.transAxes, fontsize=10,
                      va="bottom", family="monospace", bbox=box)
    txt_p = ax_p.text(0.03, 0.03, "", transform=ax_p.transAxes, fontsize=10,
                      va="bottom", family="monospace", bbox=box)

    ax_s.legend(loc="upper left")
    ax_p.legend(loc="upper left")
    fig.tight_layout()

    # ------------------------------------------------------------------
    # Animation
    # ------------------------------------------------------------------
    STEP = 2        # ft per frame
    HOLD = 20       # extra frames to hold the final picture
    frames = list(range(0, len(depths), STEP)) + [len(depths) - 1] * HOLD

    def update(k):
        d = depths[k]
        kk = k + 1                       # number of points drawn so far
        ko = min(kk, n_oil)
        kg = min(kk, n_gas)

        # curves
        l_sw.set_data(Sw_[:kk], depths[:kk])
        l_so.set_data(So_[:kk], depths[:kk])
        l_sg.set_data(Sg_[:kk], depths[:kk])
        l_pw.set_data(Pw_[:kk], depths[:kk])
        l_po.set_data(Po_[:ko], depths[:ko])
        l_pg.set_data(Pg_[:kg], depths[:kg])

        # markers
        d_sw.set_data([Sw_[k]], [d])
        d_so.set_data([So_[k]], [d])
        d_sg.set_data([Sg_[k]], [d])
        d_pw.set_data([Pw_[k]], [d])
        d_po.set_data([Po_[k]] if k < n_oil else [], [d] if k < n_oil else [])
        d_pg.set_data([Pg_[k]] if k < n_gas else [], [d] if k < n_gas else [])

        # saturation value labels follow the dots (flip to the left near the right edge)
        for lab, arr, name in sat_labels:
            lab.xy = (arr[k], d)
            lab.set_text(f"{name} = {arr[k]:.3f}")
            if arr[k] > 0.75:
                lab.set_ha("right")
                lab.xyann = (-10, lab.xyann[1])
            else:
                lab.set_ha("left")
                lab.xyann = (10, lab.xyann[1])

        # cursors
        cur_s.set_ydata([d, d])
        cur_p.set_ydata([d, d])

        # read-outs
        txt_s.set_text(f"Depth = {d:.0f} ft\n"
                       f"Sw = {Sw_[k]:.3f}\n"
                       f"So = {So_[k]:.3f}\n"
                       f"Sg = {Sg_[k]:.3f}")
        po_txt = f"{Po_[k]:.1f}" if k < n_oil else "  - "
        pg_txt = f"{Pg_[k]:.1f}" if k < n_gas else "  - "
        txt_p.set_text(f"Depth = {d:.0f} ft\n"
                       f"Pw = {Pw_[k]:.1f} psi\n"
                       f"Po = {po_txt} psi\n"
                       f"Pg = {pg_txt} psi")
        return (l_sw, l_so, l_sg, l_pw, l_po, l_pg,
                d_sw, d_so, d_sg, d_pw, d_po, d_pg,
                cur_s, cur_p, txt_s, txt_p, lab_sw, lab_so, lab_sg)

    ani = FuncAnimation(fig, update, frames=frames, interval=60, blit=False)

    ani.save("saturation_pressure_vs_depth.gif", writer=PillowWriter(fps=15), dpi=80)
    print("Saved saturation_pressure_vs_depth.gif")
    plt.show()
