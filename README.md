# Operational Range Calculator

This web application calculates and visualizes the operational ranges for three different aerosol measurement instruments:
- Differential Mobility Analyzer (DMA)
- Aerodynamic Aerosol Classifier (AAC)
- Centrifugal Particle Mass Analyzer (CPMA)

## How to Use

1. For each instrument section:
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
- Mass Resolution (R_m)

## Dependencies

All dependencies are loaded from CDN:
- Bootstrap 5.3.0
- Plotly.js 2.24.1
- Pyodide 0.23.4

No local installation is required. 
