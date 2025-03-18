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
        console.log('Starting calculation...');
        
        await initializePyodide();
        
        const Q_a = parseFloat(document.getElementById('qa').value) / 60000;
        const Q_sh = parseFloat(document.getElementById('qsh').value) / 60000;
        
        console.log('Q_a:', Q_a, 'Q_sh:', Q_sh);
        
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
                title: 'Mobility diameter, d_m (nm)',
                type: 'log',
                showgrid: true,
                gridwidth: 1,
                range: [1, 3],  // log10 range: 10^1 to 10^3 nm
                dtick: 1,       // One tick per power of 10
                ticktext: ['10', '100', '1000'],
                tickvals: [1, 2, 3]
            },
            yaxis: {
                title: 'R_B',
                type: 'log',
                showgrid: true,
                gridwidth: 1,
                range: [0, 2],  // log10 range: 10^0 to 10^2
                dtick: 1,       // One tick per power of 10
                ticktext: ['1', '10', '100'],
                tickvals: [0, 1, 2]
            },
            showlegend: true,
            legend: {
                x: 0.7,
                y: 0.9,
                xanchor: 'left',
                yanchor: 'top'
            },
            annotations: [
                {
                    x: d_i * 1e9,
                    y: Q_sh / Q_a,
                    text: 'd_m,i = ' + formatScientific(d_i * 1e9) + ' nm',
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
                    text: 'd_m,o = ' + formatScientific(d_o * 1e9) + ' nm',
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
            ],
            width: 800,
            height: 600,
            margin: {
                l: 80,
                r: 50,
                t: 50,
                b: 80
            }
        };
        
        const config = {
            displayModeBar: true,
            responsive: true
        };

        console.log('Plotting...');
        Plotly.newPlot('plot', traces, layout, config);
        console.log('Plot complete');
    } catch (error) {
        console.error('Detailed error:', error);
        alert('An error occurred: ' + error.message);
    }
}

async function calculateAAC() {
    try {
        console.log('Starting AAC calculation...');
        
        // Constants
        const P = 101325;
        const T = 298.15;
        const mu = 1.81809e-5 * Math.pow(T / 293.15, 1.5) * (293.15 + 110.4) / (T + 110.4);
        const Q_a = parseFloat(document.getElementById('aac-qa').value) / 60000;
        const Q_sh = parseFloat(document.getElementById('aac-qsh').value) / 60000;
        
        console.log('Q_a:', Q_a, 'Q_sh:', Q_sh);
        
        const Q_sh_lb = 2 / 60000;
        const Q_sh_ub = 15 / 60000;
        const Q_sh_RB = 10 / 60000;
        const r_1 = 56e-3;
        const r_2 = 60e-3;
        const L = 0.206;
        const w_lb_i = 2 * Math.PI / 60 * 200;
        const w_ub_i = 2 * Math.PI / 60 * 7000;

        const points = 100;
        const Q_sh_spa = Array.from({length: points}, (_, i) => 
            Q_sh_lb + (Q_sh_ub - Q_sh_lb) * i / (points - 1));
        const R_t = Q_sh_spa.map(q => q / Q_a);
        
        console.log('Calculating w_up...');
        const w_up = Q_sh_spa.map(Q => {
            if (Q < Q_sh_RB) {
                return Math.min(w_ub_i, 723.7 - 9.87 * 60000 * Q);
            } else {
                return Math.min(w_ub_i, 875 - 25 * 60000 * Q);
            }
        });

        const w_low = Array(points).fill(w_lb_i);

        const factor1 = w_low.map(w => (36 * mu) / (Math.PI * 1000 * Math.pow(r_1 + r_2, 2) * L * Math.pow(w, 2)));
        const factor2 = w_up.map(w => (36 * mu) / (Math.PI * 1000 * Math.pow(r_1 + r_2, 2) * L * Math.pow(w, 2)));

        console.log('Calculating d_min and d_max...');
        const d_min = Q_sh_spa.map((Q_sh_val, i) => {
            const f = d => Math.pow(d, 2) * Cc(d, P, T) - (factor2[i] * Q_sh_val);
            return bisectionMethod(f, 1e-9, 1e-7);
        });

        const d_max = Q_sh_spa.map((Q_sh_val, i) => {
            const f = d => Math.pow(d, 2) * Cc(d, P, T) - (factor1[i] * Q_sh_val);
            return bisectionMethod(f, 1e-7, 1e-5);
        });

        console.log('Plotting...');
        const trace = {
            x: [d_min, d_max],
            y: [R_t, R_t],
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
                title: {
                    text: 'Aerodynamic diameter, d_a (nm)',
                    font: { size: 14 }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1,
                range: [1, 3],  // 10 to 1000 nm
                dtick: 1
            },
            yaxis: {
                title: {
                    text: 'R_τ',
                    font: { size: 14 }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1,
                range: [0, 2],  // 1 to 100
                dtick: 1
            },
            showlegend: true,
            legend: {
                x: 0.7,
                y: 0.9,
                xanchor: 'left',
                yanchor: 'top'
            },
            width: 800,
            height: 600,
            margin: {
                l: 80,
                r: 50,
                t: 50,
                b: 80
            }
        };
        
        Plotly.newPlot('aac-plot', [trace], layout);
    } catch (error) {
        console.error('Error in calculateAAC:', error);
        alert('An error occurred during calculation. Please check the console for details.');
    }
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

// Modify the bisection method to handle the AAC equations better
function bisectionMethod(func, a, b, tolerance = 1e-10, maxIterations = 100) {
    let fa = func(a);
    let fb = func(b);
    
    if (fa * fb > 0) {
        throw new Error('Initial points do not bracket the root');
    }
    
    for (let i = 0; i < maxIterations; i++) {
        let c = (a + b) / 2;
        let fc = func(c);
        
        if (Math.abs(fc) < tolerance) {
            return c;
        }
        
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    
    return (a + b) / 2;
}

// Update the calculation of factors
const factor1 = (36 * mu) / (Math.PI * 1000 * Math.pow(r_1 + r_2, 2) * L * Math.pow(w_lb_i, 2));
const factor2 = (36 * mu) / (Math.PI * 1000 * Math.pow(r_1 + r_2, 2) * L * Math.pow(w_ub_i, 2));

// Modify the root finding functions
const d_min = Q_sh_spa.map(Q_sh_val => {
    const f = d => Math.pow(d, 2) * Cc(d, P, T) - (factor2 * Q_sh_val);
    try {
        return bisectionMethod(f, 1e-9, 1e-7);
    } catch (error) {
        console.error('Error finding d_min:', error);
        return null;
    }
});

const d_max = Q_sh_spa.map(Q_sh_val => {
    const f = d => Math.pow(d, 2) * Cc(d, P, T) - (factor1 * Q_sh_val);
    try {
        return bisectionMethod(f, 1e-7, 1e-5);
    } catch (error) {
        console.error('Error finding d_max:', error);
        return null;
    }
});

// Add error checking before plotting
if (d_min.includes(null) || d_max.includes(null)) {
    throw new Error('Failed to calculate some diameter values');
} 