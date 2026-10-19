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

![Uploading saturation_pressure_vs_depth.gif…]()



