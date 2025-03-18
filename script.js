let pyodide = null;

async function initializePyodide() {
    if (!pyodide) {
        pyodide = await loadPyodide();
        await pyodide.loadPackage(['numpy', 'scipy']);
        
        // Load the Python functions
        await pyodide.runPythonAsync(`
            import numpy as np
            from scipy.optimize import fsolve

            def Cc(d, P, T):
                alpha = 1.165 * 2
                beta = 0.483 * 2
                gamma = 0.997 / 2
                la = 67.30e-9
                lap = la * (T / 296.15) ** 2 * (101325/P) * ((110.4 + 296.15) / (T+110.4))
                return 1 + lap / d * (alpha + beta * np.exp(-gamma * d / lap))
        `);
    }
}

async function calculateDMA() {
    try {
        await initializePyodide();
        
        const Q_a = parseFloat(document.getElementById('qa').value) / 60000;
        
        const result = await pyodide.runPythonAsync(`
            def calculate_DMA(Q_a):
                P = 101325
                T = 298.15
                mu = 1.81809e-5 * (T / 293.15) ** 1.5 * (293.15 + 110.4) / (T + 110.4)
                Q_a = Q_a / 60000

                Q_sh_lb = 2 / 60000
                Q_sh_ub = 30 / 60000
                r_1 = 9.37e-3
                r_2 = 19.61e-3
                L = 0.44369
                e = 1.6e-19
                V_min = 10
                V_max = 10000

                log_r_ratio = np.log(r_2 / r_1)
                factor1 = (2 * V_min * L * e) / (3 * mu * log_r_ratio)
                factor2 = (2 * V_max * L * e) / (3 * mu * log_r_ratio)

                Q_sh_spa = np.linspace(Q_sh_lb, Q_sh_ub, 100)  # Reduced points for faster calculation
                R_B = Q_sh_spa / Q_a
                R_B_lb = Q_sh_lb / Q_a
                R_B_up = Q_sh_ub / Q_a

                d_min = np.zeros_like(Q_sh_spa)
                d_max = np.zeros_like(Q_sh_spa)

                for i, Q_sh_val in enumerate(Q_sh_spa):
                    d_min[i] = float(fsolve(lambda d: d - (factor1 / Q_sh_val) * Cc(d, P, T), 1e-8)[0])
                    d_max[i] = float(fsolve(lambda d: d - (factor2 / Q_sh_val) * Cc(d, P, T), 1e-6)[0])

                return {
                    'd_min': d_min.tolist(),
                    'd_max': d_max.tolist(),
                    'R_B': R_B.tolist(),
                    'R_B_lb': float(R_B_lb),
                    'R_B_up': float(R_B_up)
                }

            calculate_DMA(${Q_a})
        `);

        console.log('Calculation result:', result);  // Debug output

        // Convert diameter values from meters to nanometers
        const d_min_nm = result.d_min.map(d => d * 1e9);
        const d_max_nm = result.d_max.map(d => d * 1e9);
        const R_B = result.R_B;
        const R_B_lb = result.R_B_lb;
        const R_B_up = result.R_B_up;

        // Format numbers to scientific notation with 2 decimal places
        function formatScientific(num) {
            return num.toExponential(2);
        }

        const traces = [
            {
                x: d_min_nm,
                y: R_B,
                mode: 'lines',
                line: { color: 'red' },
                showlegend: false
            },
            {
                x: d_max_nm,
                y: R_B,
                mode: 'lines',
                line: { color: 'red' },
                showlegend: false
            },
            {
                x: [d_min_nm[0], d_max_nm[0]],
                y: [R_B_lb, R_B_lb],
                mode: 'lines',
                line: { color: 'red' },
                showlegend: false
            },
            {
                x: [d_min_nm[d_min_nm.length-1], d_max_nm[d_max_nm.length-1]],
                y: [R_B_up, R_B_up],
                mode: 'lines',
                line: { color: 'red' },
                showlegend: false
            },
            {
                x: [d_i * 1e9, d_o * 1e9],
                y: [Q_sh / Q_a, Q_sh / Q_a],
                mode: 'lines',
                name: 'Selected Q_sh',
                line: { color: 'blue' }
            }
        ];
        
        const layout = {
            title: 'DMA Operational Range',
            xaxis: {
                title: {
                    text: 'Mobility diameter, $d_{\\mathrm{m}}$ (nm)',
                    font: {
                        size: 14
                    }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1
            },
            yaxis: {
                title: {
                    text: '$R_{\\mathrm{B}}$',
                    font: {
                        size: 14
                    }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1
            },
            showlegend: true,
            annotations: [
                {
                    x: d_i * 1e9,
                    y: Q_sh / Q_a,
                    text: `$d_{\\mathrm{m,l}}$ = ${formatScientific(d_i * 1e9)} nm`,
                    showarrow: true,
                    arrowhead: 2,
                    arrowsize: 1,
                    arrowwidth: 1,
                    ax: -40,
                    ay: -40,
                    font: {
                        color: 'blue'
                    }
                },
                {
                    x: d_o * 1e9,
                    y: Q_sh / Q_a,
                    text: `$d_{\\mathrm{m,u}}$ = ${formatScientific(d_o * 1e9)} nm`,
                    showarrow: true,
                    arrowhead: 2,
                    arrowsize: 1,
                    arrowwidth: 1,
                    ax: 40,
                    ay: -40,
                    font: {
                        color: 'blue'
                    }
                }
            ]
        };
        
        // Add MathJax config to enable LaTeX rendering
        const config = {
            mathjax: 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_SVG'
        };

        Plotly.newPlot('plot', traces, layout, config);
    } catch (error) {
        console.error('Error in calculateDMA:', error);
        alert('An error occurred during calculation. Please check the console for details.');
    }
}

async function calculateAAC() {
    await initializePyodide();
    
    const Q_a = parseFloat(document.getElementById('aac-qa').value);
    const Q_sh = parseFloat(document.getElementById('aac-qsh').value);
    
    const result = await pyodide.runPythonAsync(`calculate_AAC(${Q_a}, ${Q_sh})`);
    
    const trace = {
        x: [result.d_min, result.d_max],
        y: [result.R_t, result.R_t],
        mode: 'lines',
        name: 'AAC Range',
        line: {
            color: 'blue',
            width: 2
        }
    };
    
    const layout = {
        title: 'AAC Operational Range',
        xaxis: {
            title: 'Particle diameter (nm)',
            type: 'log'
        },
        yaxis: {
            title: 'R_t',
            type: 'log'
        }
    };
    
    Plotly.newPlot('aac-plot', [trace], layout);
}

async function calculateCPMA() {
    await initializePyodide();
    
    const Q_a = parseFloat(document.getElementById('cpma-qa').value);
    
    const result = await pyodide.runPythonAsync(`calculate_CPMA(${Q_a})`);
    
    const trace = {
        x: [result.d_min, result.d_max],
        y: [result.R_m, result.R_m],
        mode: 'lines',
        name: 'CPMA Range',
        line: {
            color: 'green',
            width: 2
        }
    };
    
    const layout = {
        title: 'CPMA Operational Range',
        xaxis: {
            title: 'Particle diameter (nm)',
            type: 'log'
        },
        yaxis: {
            title: 'R_m',
            type: 'log'
        }
    };
    
    Plotly.newPlot('cpma-plot', [trace], layout);
}

// Calculate initial plot on page load
document.addEventListener('DOMContentLoaded', calculateDMA); 