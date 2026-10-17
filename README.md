# Brooks-Corey Saturation Model

This project implements a simple Brooks-Corey reservoir saturation model and plots how water, oil, and gas saturations vary with depth.

## What the model does

The script estimates phase behavior across a vertical depth interval using the following inputs:

- `WOC`: water-oil contact depth
- `GOC`: gas-oil contact depth
- fluid densities for the aqueous and oil phases
- capillary entry pressure and fluid property parameters

For each depth, the model calculates:

- phase pressures for water, oil, and gas
- saturation values for water, oil, and gas
- vertical trends from 7300 ft to 7550 ft

The result is a depth profile showing how fluid saturations shift around the contact levels.

<img width="640" height="480" alt="brooks_corey_animation" src="https://github.com/user-attachments/assets/886c008b-216e-4f16-8090-62c700015e67" />


