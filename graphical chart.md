# Graphical Parameters for Python Plots 📊

This document outlines the standard graphical parameters to use for creating plots in Python for our ChE-320 projects.

1. [Line Plot Parameters](#line-plot-parameters-)
2. [Scatter Plot Parameters](#scatter-plot-parameters-)
3. [Grids Parameters](#grids-parameters-)
4. [Color Scheme](#color-scheme-)
5. [Graph Size and Layout](#graph-size-and-layout-)
6. [Font and Text Parameters](#font-and-text-parameters-)
7. [Plot annotations](#plot-annotations-%EF%B8%8F)

## Line Plot Parameters 📈

- **Line Width**: 
  - Use a line width `linewidth = '1.0'`.

- **Line Style**:
  - Use solid lines `linestyle = '-'` for data line representing a model.
  - Use dashed lines `linestyle = '--'` for data line representing samples.
  - Use dotted lines `linestyle = ':'` for data line highlighting important values.

## Scatter Plot Parameters 🔍

- **Marker Size**:
  - Use marker size of `s = 10.0`.

- **Marker Style**:
  - Use circular markers `marker = 'o'`.
 
## Grids Parameters 📐

- **Line Width**: 
  - Use a line width of `0.25`.

- **Line Style**:
  - Use solid lines `linestyle = '-'`.

## Color Scheme 🎨

- **1 Color**:
  - Use `color = 'black'` if no other same linestyle are used on the same plot.

- **2 Colors**:
  - Use `color = 'red'` and `color = 'blue'`.

## Graph Size and Layout 📋

- **Figure Size**:
  - Use `plt.rcParams['figure.figsize'] = [8, 3]` at the beginning of your code.
  - For figure with multiple plots, use `fig, axs = plt.subplots(2, 2, figsize=(8, nbl*3))` where nbl is the number of subplots lines of the figure.
 
- **Figure Layout**:
  - Use `plt.rcParams['figure.constrained_layout.use'] = True` at the beginning of your code.
 
- **Axes Limits and Ticks**:

## Font and Text Parameters 📝

- **Font**:
  - Use `not defined`.

- **Font Size**:
  - Use `plt.rcParams.update({'font.size': 10})` at the beginning of your code.
 
## Plot annotations 🖋️

- **Title**:

- **Axes**:

- **Legend**:

- **Annotation**:

- **Units**: