# Graphical Parameters for Python Plots 📊

This document outlines the standard graphical parameters to use for creating plots in Python for our ChE-320 projects.

## Line Plot Parameters 📈

- **Line Width**: 
  - Use a line width of `1.0`.

- **Line Style**:
  - Use solid lines (`linestyle = '-'`) for data line representing a model.
  - Use dashed lines (`linestyle = '--'`) for data line representing samples.

## Scatter Plot Parameters 🔍

- **Marker Size**:
  - Use a marker size of `10.0`.

- **Marker Style**:
  - Use circular markers (`o`).

## Color Scheme 🎨

- **1 Color**:
  - Use `black` if no other same linestyle are used on the same plot.

- **2 Colors**:
  - Use `red` and `blue`.

## Graph Size and Layout 📐

- **Figure Size**:
  - Use `plt.rcParams['figure.figsize'] = [8, 3]` (width x height) at the beginning of your code.
  - For figure with multiple plots, use `fig, axs = plt.subplots(2, 2, figsize=(8, nbl*3))` where nbl is the number of subplots lines of the figure.

## Font and Text Parameters 📝

- **Font Size**:
  - Use `plt.rcParams.update({'font.size': 10})` at the beginning of your code.
 
## Units 📏

