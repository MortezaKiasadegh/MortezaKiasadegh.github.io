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

        // Format numbers to scientific notation with 2 decimal places
        function formatScientific(num) {
            return num.toExponential(2);
        }

        const traces = [
            {
                x: d_min_nm,
                y: result.R_B,
                mode: 'lines',
                name: 'Lower bound',
                line: { color: 'red' }
            },
            {
                x: d_max_nm,
                y: result.R_B,
                mode: 'lines',
                name: 'Upper bound',
                line: { color: 'red' }
            },
            {
                x: [d_min_nm[0], d_max_nm[0]],
                y: [result.R_B_lb, result.R_B_lb],
                mode: 'lines',
                name: 'Lower R_B',
                line: { color: 'red' }
            },
            {
                x: [d_min_nm[d_min_nm.length-1], d_max_nm[d_max_nm.length-1]],
                y: [result.R_B_up, result.R_B_up],
                mode: 'lines',
                name: 'Upper R_B',
                line: { color: 'red' }
            },
            {
                x: [d_i * 1e9, d_o * 1e9],
                y: [Q_sh / Q_a, Q_sh / Q_a],
                mode: 'lines+text',
                name: 'Selected Q_sh',
                line: { color: 'blue' },
                text: [`d_i = ${formatScientific(d_i * 1e9)} nm`, `d_o = ${formatScientific(d_o * 1e9)} nm`],
                textposition: ['bottom left', 'bottom right'],
                textfont: {
                    family: 'Arial',
                    size: 12,
                    color: 'blue'
                }
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
    try {
        // Constants
        const P = 101325;
        const T = 298.15;
        const mu = 1.81809e-5 * Math.pow(T / 293.15, 1.5) * (293.15 + 110.4) / (T + 110.4);
        const Q_a = parseFloat(document.getElementById('aac-qa').value) / 60000;
        const Q_sh = parseFloat(document.getElementById('aac-qsh').value) / 60000;
        
        // Classifier properties
        const r_1 = 56e-3;
        const r_2 = 60e-3;
        const L = 0.206;
        const Q_sh_lb = 2 / 60000;
        const Q_sh_ub = 15 / 60000;
        const Q_sh_RB = 10 / 60000;
        const w_lb_i = 2 * Math.PI / 60 * 200;
        const w_ub_i = 2 * Math.PI / 60 * 7000;

        console.log('Starting AAC calculation...'); // Debug log

        // Create arrays for Q_sh values using logarithmic spacing
        const points = 100;
        const Q_sh_spa = [];
        for (let i = 0; i < points; i++) {
            const t = i / (points - 1);
            Q_sh_spa.push(Q_sh_lb * Math.pow(Q_sh_ub/Q_sh_lb, t));
        }

        const R_t = Q_sh_spa.map(q => q / Q_a);
        const R_t_lb = Q_sh_lb / Q_a;
        const R_t_up = Q_sh_ub / Q_a;

        // Calculate w_up for each Q_sh value
        const w_low = Array(points).fill(w_lb_i);
        const w_up = Q_sh_spa.map(Q => {
            if (Q < Q_sh_RB) {
                return Math.min(w_ub_i, 723.7 - 9.87 * 60000 * Q);
            } else {
                return Math.min(w_ub_i, 875 - 25 * 60000 * Q);
            }
        });

        // Calculate factors
        const factor1 = w_low.map(w => (36 * mu) / (Math.PI * 1000 * Math.pow(r_1 + r_2, 2) * L * Math.pow(w, 2)));
        const factor2 = w_up.map(w => (36 * mu) / (Math.PI * 1000 * Math.pow(r_1 + r_2, 2) * L * Math.pow(w, 2)));

        console.log('Calculating d_min and d_max...'); // Debug log

        // Calculate d_min and d_max with better initial guesses
        const d_min = Q_sh_spa.map((Q_sh_val, i) => {
            const f = d => Math.pow(d, 2) * Cc(d, P, T) - (factor2[i] * Q_sh_val);
            // Use wider initial bracket for d_min
            return bisectionMethod(f, 1e-10, 1e-6);
        });

        const d_max = Q_sh_spa.map((Q_sh_val, i) => {
            const f = d => Math.pow(d, 2) * Cc(d, P, T) - (factor1[i] * Q_sh_val);
            // Use wider initial bracket for d_max
            return bisectionMethod(f, 1e-8, 1e-4);
        });

        // Calculate specific points for input Q_sh with better initial guesses
        const idx = Q_sh_spa.findIndex(q => Math.abs(q - Q_sh) < Q_sh * 0.01);
        const d_i = bisectionMethod(
            d => Math.pow(d, 2) * Cc(d, P, T) - (factor2[idx] * Q_sh),
            1e-10, 1e-6
        );
        const d_o = bisectionMethod(
            d => Math.pow(d, 2) * Cc(d, P, T) - (factor1[idx] * Q_sh),
            1e-8, 1e-4
        );

        console.log('Creating plot...'); // Debug log

        // Create plot
        const traces = [
            {
                x: d_min.map(d => d * 1e9),
                y: R_t,
                mode: 'lines',
                name: 'Lower bound',
                line: { color: 'red' }
            },
            {
                x: d_max.map(d => d * 1e9),
                y: R_t,
                mode: 'lines',
                name: 'Upper bound',
                line: { color: 'red' }
            },
            {
                x: [d_min[0] * 1e9, d_max[0] * 1e9],
                y: [R_t_lb, R_t_lb],
                mode: 'lines',
                name: 'Lower R_τ',
                line: { color: 'red' }
            },
            {
                x: [d_min[d_min.length-1] * 1e9, d_max[d_max.length-1] * 1e9],
                y: [R_t_up, R_t_up],
                mode: 'lines',
                name: 'Upper R_τ',
                line: { color: 'red' }
            },
            {
                x: [d_i * 1e9, d_o * 1e9],
                y: [Q_sh / Q_a, Q_sh / Q_a],
                mode: 'lines',
                name: 'Selected Q_sh',
                line: { color: 'green' }
            }
        ];

        const layout = {
            title: 'AAC Operational Range',
            xaxis: {
                title: {
                    text: 'Aerodynamic diameter, $d_{\\mathrm{a}}$ (nm)',
                    font: { size: 14 }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1
            },
            yaxis: {
                title: {
                    text: '$R_{\\tau}$',
                    font: { size: 14 }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1
            },
            showlegend: true,
            grid: { pattern: 'independent' }
        };

        console.log('Plotting...'); // Debug log

        Plotly.newPlot('aac-plot', traces, layout);

    } catch (error) {
        console.error('Detailed error in calculateAAC:', error);
        console.log('Input values:', {
            Q_a: document.getElementById('aac-qa').value,
            Q_sh: document.getElementById('aac-qsh').value
        });
        alert('An error occurred during AAC calculation. Please check the console for details.');
    }
}

async function calculateCPMA() {
    try {
        // Constants
        const P = 101325;
        const T = 298.15;
        const mu = 1.81809e-5 * Math.pow(T / 293.15, 1.5) * (293.15 + 110.4) / (T + 110.4);
        const Q_a = parseFloat(document.getElementById('cpma-qa').value) / 60000;
        const R_m_inp = parseFloat(document.getElementById('cpma-rm').value);

        // Classifier parameters
        const r_1 = 60e-3;
        const r_2 = 61e-3;
        const r_c = 60.5e-3;
        const L = 0.2;
        const e = 1.6e-19;
        const V_min = 0.1;
        const V_max = 1000;
        const w_lb_i = 2 * Math.PI / 60 * 200;
        const w_ub_i = 2 * Math.PI / 60 * 12000;
        const k = 1000 * Math.PI / 6;
        const D_m = 3;

        // Create R_m array
        const points = 100;
        const R_m_spa = Array.from({length: points}, (_, i) => {
            const t = i / (points - 1);
            return Math.pow(10, Math.log10(0.001) * (1-t) + Math.log10(200) * t);
        });

        // Precompute factors
        const log_r_ratio = Math.log(r_2 / r_1);
        const factor1 = e * V_min / (k * Math.pow(r_c, 2) * log_r_ratio);
        const factor2 = e * V_max / (k * Math.pow(r_c, 2) * log_r_ratio);
        const factor3 = 3 * mu * Q_a / (2 * k * Math.pow(r_c, 2) * L);

        // Residual function for optimization
        function residual(d, R_m_val, factor_v, is_min) {
            const d_m_max = Math.pow((R_m_val + 1) / R_m_val, 1 / D_m) * d;
            const factor = factor_v === 'min' ? factor1 : factor2;
            
            const w_guess = Math.sqrt(factor / Math.pow(d, D_m));
            const w = Math.min(Math.max(w_guess, w_lb_i), w_ub_i);
            
            return (Math.pow(d, D_m) - Math.pow(d_m_max, D_m) + 
                   (factor3 / Math.pow(w, 2)) * (d_m_max / Cc(d_m_max, P, T))) / 
                   Math.abs(Math.pow(d, D_m));
        }

        // Function to find root using bisection
        function findRoot(R_m_val, factor_v, is_min, lower, upper) {
            const f = d => residual(d, R_m_val, factor_v, is_min);
            return bisectionMethod(f, lower, upper);
        }

        // Calculate boundaries for input R_m
        const d_i_1 = findRoot(R_m_inp, 'min', true, 1e-9, 1e-7);
        const d_o_1 = findRoot(R_m_inp, 'min', false, 1e-7, 1e-5);
        const d_i_2 = findRoot(R_m_inp, 'max', true, 1e-9, 1e-7);
        const d_o_2 = findRoot(R_m_inp, 'max', false, 1e-7, 1e-5);

        const d_i = Math.max(d_i_1, d_i_2);
        const d_o = Math.min(d_o_1, d_o_2);

        // Calculate boundaries for R_m sweep
        const d_min_1 = R_m_spa.map(R_m => findRoot(R_m, 'min', true, 1e-9, 1e-7));
        const d_max_1 = R_m_spa.map(R_m => findRoot(R_m, 'min', false, 1e-7, 1e-5));
        const d_min_2 = R_m_spa.map(R_m => findRoot(R_m, 'max', true, 1e-9, 1e-7));
        const d_max_2 = R_m_spa.map(R_m => findRoot(R_m, 'max', false, 1e-7, 1e-5));

        // Combine results
        const d_min = d_min_1.map((d, i) => Math.max(d, d_min_2[i]));
        const d_max = d_max_1.map((d, i) => Math.min(d, d_max_2[i]));

        // Filter valid points
        const valid_points = d_min.map((d, i) => ({
            d_min: d,
            d_max: d_max[i],
            R_m: R_m_spa[i]
        })).filter(p => Math.abs(p.d_min - p.d_max) > 1e-8);

        // Create plot
        const traces = [
            {
                x: valid_points.map(p => p.d_min * 1e9),
                y: valid_points.map(p => p.R_m),
                mode: 'lines',
                name: 'Lower bound',
                line: { color: 'red' }
            },
            {
                x: valid_points.map(p => p.d_max * 1e9),
                y: valid_points.map(p => p.R_m),
                mode: 'lines',
                name: 'Upper bound',
                line: { color: 'red' }
            },
            {
                x: [d_i * 1e9, d_o * 1e9],
                y: [R_m_inp, R_m_inp],
                mode: 'lines',
                name: 'Selected R_m',
                line: { color: 'green' }
            }
        ];

        const layout = {
            title: 'CPMA Operational Range',
            xaxis: {
                title: {
                    text: 'Mobility diameter, $d_{\\mathrm{m}}$ (nm)',
                    font: { size: 14 }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1
            },
            yaxis: {
                title: {
                    text: '$R_{\\mathrm{m}}$',
                    font: { size: 14 }
                },
                type: 'log',
                showgrid: true,
                gridwidth: 1
            },
            showlegend: true,
            grid: { pattern: 'independent' }
        };

        const config = {
            mathjax: 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_SVG'
        };

        Plotly.newPlot('cpma-plot', traces, layout, config);

    } catch (error) {
        console.error('Error in calculateCPMA:', error);
        alert('An error occurred during CPMA calculation. Please check the console for details.');
    }
}

function bisectionMethod(func, a, b, tolerance = 1e-10, maxIterations = 100) {
    try {
        let left = a;
        let right = b;
        let fa = func(left);
        let fb = func(right);
        
        // Check if initial points bracket a root
        if (fa * fb >= 0) {
            console.log('Initial points:', left, right);
            console.log('Function values:', fa, fb);
            throw new Error('Initial points do not bracket the root');
        }

        for (let i = 0; i < maxIterations; i++) {
            let mid = (left + right) / 2;
            let fmid = func(mid);
            
            if (Math.abs(fmid) < tolerance) {
                return mid;
            }
            
            if (fa * fmid < 0) {
                right = mid;
                fb = fmid;
            } else {
                left = mid;
                fa = fmid;
            }
        }
        
        return (left + right) / 2;
    } catch (error) {
        console.error('Error in bisectionMethod:', error);
        throw error;
    }
}

// Calculate initial plot on page load
document.addEventListener('DOMContentLoaded', () => {
    calculateDMA();
    calculateAAC();
    calculateCPMA();
}); 