This project implements a simple Brooks-Corey reservoir saturation model and plots how water, oil, and gas saturations vary with depth.

## What the model does

The script estimates phase behavior across a vertical depth interval using the following inputs:

- `WOC`: water-oil contact depth
- `GOC`: gas-oil contact depth
- `FWL` - free water level
- `FGL` - free gas level
- fluid densities for the aqueous and oil phases
- capillary entry pressure and fluid property parameters

For each depth, the model calculates:

- phase pressures for water, oil, and gas
- saturation values for water, oil, and gas
- vertical trends from 7300 ft to 7550 ft

The result is a depth profile showing how fluid saturations shift around the contact levels.

<img width="960" height="640" alt="saturation_pressure_vs_depth" src="https://github.com/user-attachments/assets/f2d4b272-4f1a-4919-89ad-f3b261929efe" />








A 3D scene:

<img width="800" height="492" alt="fwl_vs_woc_animation" src="https://github.com/user-attachments/assets/722aa914-ce02-4521-bb3f-055cc7bfe781" />


---------------------------------------------------------------------------------------------------------------------------------------------------------------------
Used formulas:

<img width="3200" height="2120" alt="image" src="https://github.com/user-attachments/assets/c4838b25-5765-4719-b60b-2b737dabcfec" />





