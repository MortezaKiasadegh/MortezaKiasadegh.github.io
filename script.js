let pyodide = null;

async function initializePyodide() {
    if (!pyodide) {
        pyodide = await loadPyodide();
        await pyodide.loadPackage(['numpy', 'scipy']);
        
        // Load the Python functions
        await pyodide.runPythonAsync(`
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import least_squares

def Cc(d, P, T):
    alpha = 1.165 * 2
    beta = 0.483 * 2
    gamma = 0.997 / 2
    la = 67.30e-9
    lap = la * (T / 296.15) ** 2 * (101325/P) * ((110.4 + 296.15) / (T+110.4))
    return 1 + lap / d * (alpha + beta * np.exp(-gamma * d / lap))

def calculate_DMA(Q_a, Q_sh):
    P = 101325
    T = 298.15
    mu = 1.81809e-5 * (T / 293.15) ** 1.5 * (293.15 + 110.4) / (T + 110.4)
    Q_a = Q_a / 60000
    Q_sh = Q_sh / 60000
    r_1 = 9.37e-3
    r_2 = 19.61e-3
    L = 0.44369
    e = 1.6e-19
    V_min = 10
    V_max = 10000

    log_r_ratio = np.log(r_2 / r_1)
    factor1 = (2 * V_min * L * e) / (3 * mu * log_r_ratio)
    factor2 = (2 * V_max * L * e) / (3 * mu * log_r_ratio)

    def f1(d_min):
        return d_min - (factor1 / Q_sh) * Cc(d_min, P, T)

    def f2(d_max):
        return d_max - (factor2/Q_sh) * Cc(d_max, P, T)

    d_min = float(fsolve(f1, 1e-8))
    d_max = float(fsolve(f2, 1e-6))
    R_B = Q_sh / Q_a

    return {
        'd_min': d_min * 1e9,
        'd_max': d_max * 1e9,
        'R_B': R_B
    }

def calculate_AAC(Q_a, Q_sh):
    P = 101325
    T = 298.15
    mu = 1.81809e-5 * (T / 293.15) ** 1.5 * (293.15 + 110.4) / (T + 110.4)
    Q_a = Q_a / 60000
    Q_sh = Q_sh / 60000
    r_1 = 56e-3
    r_2 = 60e-3
    L = 0.206
    w_lb_i = 2 * np.pi / 60 * 200
    w_ub_i = 2 * np.pi / 60 * 7000

    R_t = Q_sh / Q_a
    
    if Q_sh < (10 / 60000):
        w_up = min(w_ub_i, 723.7 - 9.87 * 60000 * Q_sh)
    else:
        w_up = min(w_ub_i, 875 - 25 * 60000 * Q_sh)
    
    w_low = w_lb_i

    factor1 = (36 * mu) / (np.pi * 1000 * (r_1 + r_2)**2 * L * (w_low**2))
    factor2 = (36 * mu) / (np.pi * 1000 * (r_1 + r_2)**2 * L * (w_up**2))

    def f1(d_min):
        return d_min**2 * Cc(d_min, P, T) - (factor2 * Q_sh)

    def f2(d_max):
        return d_max**2 * Cc(d_max, P, T) - (factor1 * Q_sh)

    d_min = float(fsolve(f1, 1e-8))
    d_max = float(fsolve(f2, 1e-6))

    return {
        'd_min': d_min * 1e9,
        'd_max': d_max * 1e9,
        'R_t': R_t
    }

def calculate_CPMA(Q_a):
    P = 101325
    T = 298.15
    mu = 1.81809e-5 * (T / 293.15) ** 1.5 * (293.15 + 110.4) / (T + 110.4)
    Q_a = Q_a / 60000
    r_1 = 60e-3
    r_2 = 61e-3
    r_c = 60.5e-3
    L = 0.2
    e = 1.6e-19
    V_min = 0.1
    V_max = 1000
    w_lb_i = 2 * np.pi / 60 * 200
    w_ub_i = 2 * np.pi / 60 * 12000
    k = 1000 * np.pi / 6
    D_m = 3

    R_m = 0.1  # Example value
    factor1 = e * V_min / (k * r_c**2 * np.log(r_2 / r_1))
    factor2 = e * V_max / (k * r_c**2 * np.log(r_2 / r_1))
    factor3 = 3 * mu * Q_a / (2 * k * r_c**2 * L)

    def f1(d_min):
        d_m_max = ((R_m + 1) / R_m)**(1 / D_m) * d_min
        w_max_i = min(w_ub_i, (factor1 / (d_min ** D_m)) ** 0.5)
        w_max = max(w_max_i, w_lb_i)
        return d_min ** D_m - d_m_max ** D_m + (factor3 / w_max ** 2) * (d_m_max / Cc(d_m_max, P, T))

    def f2(d_max):
        d_m_max = ((R_m + 1) / R_m)**(1 / D_m) * d_max
        w_min_i = max(w_lb_i, (factor1 / (d_max ** D_m)) ** 0.5)
        w_min = min(w_min_i, w_ub_i)
        return d_max ** D_m - d_m_max ** D_m + (factor3 / w_min ** 2) * (d_m_max / Cc(d_m_max, P, T))

    d_min = float(fsolve(f1, 1e-8))
    d_max = float(fsolve(f2, 1e-6))

    return {
        'd_min': d_min * 1e9,
        'd_max': d_max * 1e9,
        'R_m': R_m
    }
        `);
    }
}

async function calculateDMA() {
    await initializePyodide();
    
    const Q_a = parseFloat(document.getElementById('dma-qa').value);
    const Q_sh = parseFloat(document.getElementById('dma-qsh').value);
    
    const result = await pyodide.runPythonAsync(`calculate_DMA(${Q_a}, ${Q_sh})`);
    
    const trace = {
        x: [result.d_min, result.d_max],
        y: [result.R_B, result.R_B],
        mode: 'lines',
        name: 'DMA Range',
        line: {
            color: 'red',
            width: 2
        }
    };
    
    const layout = {
        title: 'DMA Operational Range',
        xaxis: {
            title: 'Particle diameter (nm)',
            type: 'log'
        },
        yaxis: {
            title: 'R_B',
            type: 'log'
        }
    };
    
    Plotly.newPlot('dma-plot', [trace], layout);
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