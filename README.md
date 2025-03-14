# Operational Range Calculator

This web application calculates and visualizes the operational ranges for three different aerosol measurement instruments:
- Differential Mobility Analyzer (DMA)
- Aerodynamic Aerosol Classifier (AAC)
- Centrifugal Particle Mass Analyzer (CPMA)

## Features

- Interactive input fields for flow rates
- Real-time calculation of operational ranges
- Logarithmic plots showing the operational ranges
- Modern and responsive user interface

## How to Use

1. Open `index.html` in a modern web browser (Chrome, Firefox, Safari, or Edge recommended)
2. For each instrument section:
   - Enter the required flow rate values
   - Click the "Calculate" button
   - View the resulting plot showing the operational range

### Input Parameters

#### DMA Section
- Aerosol Flow Rate (Q_a) [L/min]
- Sheath Flow Rate (Q_sh) [L/min]

#### AAC Section
- Aerosol Flow Rate (Q_a) [L/min]
- Sheath Flow Rate (Q_sh) [L/min]

#### CPMA Section
- Aerosol Flow Rate (Q_a) [L/min]

## Technical Details

The application uses:
- Pyodide for running Python code in the browser
- NumPy and SciPy for calculations
- Plotly.js for interactive plotting
- Bootstrap for responsive layout

## Dependencies

All dependencies are loaded from CDN:
- Bootstrap 5.3.0
- Plotly.js 2.24.1
- Pyodide 0.23.4

No local installation is required. 