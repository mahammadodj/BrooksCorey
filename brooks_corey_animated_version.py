import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import text
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from pathlib import Path

class ReservoirModel:
    def __init__(self, ro_aqu, ro_oleic, WOC, GOC, M_gas=19, Z=1.2, R=10.73, T=175, sigma_ow=20, sigma_og=30, sigma_gw=50, g=1/144, Swirr=0.2, _lambda=2, Pe_ow=3.5):
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

    def calculate_p_contacts(self):
        Pw_WOC = 14.7 + self.ro_aqu * self.g * self.WOC
        Po_WOC = Pw_WOC + self.Pe_ow

        Pw_GOC = Pw_WOC + self.ro_aqu * self.g * (self.GOC - self.WOC)
        Po_GOC = Po_WOC + self.ro_oleic * self.g * (self.GOC - self.WOC)

        P_rep_gas = ((Po_GOC + 99) // 100) * 100
        ro_g = (P_rep_gas * self.M_gas) / (self.Z * self.R * (self.T + 460))

        Pe_og = self.sigma_og * self.Pe_ow / self.sigma_ow
        Pg_GOC = (Po_GOC + Pe_og)

        return Pw_WOC, Po_WOC, Pw_GOC, Po_GOC, Pg_GOC, ro_g

    def calculate_pressures(self, Depth):
        Pw_WOC, Po_WOC, Pw_GOC, Po_GOC, Pg_GOC, ro_g = self.calculate_p_contacts()
        Po = Po_WOC + self.ro_oleic * self.g * (Depth - self.WOC)
        Pw = Pw_WOC + self.ro_aqu * self.g * (Depth - self.WOC)
        Pg = Pg_GOC + ro_g * self.g * (Depth - self.GOC)
        return Po, Pw, Pg
    
    def calculate_entry_pressures(self):
        Pe_ow = self.Pe_ow
        Pe_og = Pe_ow*self.sigma_og/self.sigma_ow
        Pe_wg = Pe_ow*self.sigma_gw/self.sigma_ow
        
        return Pe_ow, Pe_og, Pe_wg
    
    def calculate_saturations(self, Depth):
        Po, Pw, Pg = self.calculate_pressures(Depth)
        Pe_ow, Pe_og, Pe_wg = self.calculate_entry_pressures()

        Sw = self.Swirr + (1 - self.Swirr) * ((Po - Pw) / self.Pe_ow) ** (-self._lambda)
        So = 1 - Sw
        Sg = 0

        if Pg - Po >= Pe_og:
            Sw = self.Swirr + (1 - self.Swirr) * ((Po - Pw) / (Pe_ow))**(-self._lambda)
            Sl = self.Swirr + (1 - self.Swirr) * ((Pg - Po) / (Pe_ow * self.sigma_og / self.sigma_ow))**(-self._lambda)
            So = Sl - Sw
            Sg = 1 - Sl
        
        elif Po-Pw >= Pe_ow:
            Sw = self.Swirr + (1 - self.Swirr)*((Po - Pw) / Pe_ow)**(-self._lambda)
            So = 1 - Sw
            Sg = 0
        else:
            Sw = 1
            So = 0
            Sg = 0
            
        # Ensure saturations are within range [0, 1]
        Sw = max(0, min(1, Sw))
        So = max(0, min(1, So))
        Sg = max(0, min(1, Sg))

        return Sw, So, Sg

def animate_saturation_plot():
    WOC = 7475
    GOC = 7450
    output_path = Path(__file__).with_name('brooks_corey_animation.gif')

    reservoir = ReservoirModel(WOC=WOC, GOC=GOC, ro_aqu=66.96, ro_oleic=53)
    depths = list(range(7300, 7550))
    Sw_ = []
    So_ = []
    Sg_ = []
    depth_ = []

    for depth in depths:
        Sw, So, Sg = reservoir.calculate_saturations(Depth=depth)
        Sw_.append(Sw)
        So_.append(So)
        Sg_.append(Sg)
        depth_.append(depth)

    fig, ax = plt.subplots()
    (water_line,) = ax.plot([], [], color='blue', label='Water Saturation')
    (oil_line,) = ax.plot([], [], color='red', label='Oil Saturation')
    (gas_line,) = ax.plot([], [], color='g', label='Gas Saturation')

    ax.set_ylim(7550, 7300)
    ax.set_xlim(0, 1)
    ax.axhline(y=GOC, linestyle='--', color='black')
    text(0.9, GOC - 6, 'GOC', verticalalignment='center', fontsize=12)
    ax.axhline(y=WOC, linestyle='--', color='black')
    text(0.9, WOC - 6, 'WOC', verticalalignment='center', fontsize=12)
    ax.set_xlabel('Depths [ft]')
    ax.set_ylabel('Phase Saturation')
    ax.set_title('Fluid Saturation vs Depth')
    ax.legend()

    def update(frame):
        end = frame + 1
        water_line.set_data(Sw_[:end], depth_[:end])
        oil_line.set_data(So_[:end], depth_[:end])
        gas_line.set_data(Sg_[:end], depth_[:end])
        return water_line, oil_line, gas_line

    anim = FuncAnimation(fig, update, frames=len(depth_), interval=25, blit=True, repeat=True)
    anim.save(output_path, writer=PillowWriter(fps=40))
    plt.show()
    return anim


if __name__ == '__main__':
    animate_saturation_plot()