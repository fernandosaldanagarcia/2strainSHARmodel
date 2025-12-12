# Two-Strain SHAR Model: Code and Analysis

This repository contains code used in the following paper:

**Saldaña, F., Stollenwerk, N., & Aguiar, M. (2024).** Modelling COVID-19 mutant dynamics: understanding the interplay between viral evolution and disease transmission dynamics. *Royal Society Open Science*, 11(10), 240919.

---

## System Requirements

- **Python**: 3.9 or higher
- **Operating System**: macOS, Linux, or Windows with Python installed
- **Jupyter**: For running notebook files

## Required Python Libraries

### Core Scientific Computing
- **numpy**: Numerical computations and array operations
- **scipy**: Scientific functions, especially `scipy.integrate.odeint` for solving ODEs
- **matplotlib**: Plotting and visualization (version 3.x)

### Specialized Analysis Libraries
- **SALib**: Sensitivity analysis library for Sobol indices computation
  - Used in `AmSobol.ipynb` for global sensitivity analysis
  - Install: `pip install SALib`

### Stochastic Modeling
- **GillesPy2**: Gillespie algorithm implementation for stochastic simulations
  - Required for `two_strainG.py` and stochastic simulations
  - Install: `pip install gillespy2`
  - Includes solvers: `TauLeapingSolver` and `ODESolver`

## Installation Instructions

### Option 1: Using Virtual Environment (Recommended)

```bash
# Create a virtual environment
python -m venv shar_env
source shar_env/bin/activate  # On Windows: shar_env\Scripts\activate

# Install required packages
pip install numpy scipy matplotlib SALib gillespy2

# Install Jupyter for notebooks
pip install jupyter
```

### Option 2: Using Conda

```bash
conda create -n shar_env python=3.9
conda activate shar_env
conda install numpy scipy matplotlib
pip install SALib gillespy2 jupyter
```

---

## File Descriptions and Usage

### 1. `equilibriaSHARimport.ipynb`

**Purpose**: Computes and plots endemic equilibria for the SHAR model with pathogen import

**Key Functions**:
- `Hstar(params)`: Computes $H^*$ (hospitalized wild-type equilibrium)
- `Astar(params)`: Computes $A^*$ (asymptomatic wild-type equilibrium)

**Input Parameters** (in order): $N, \gamma, \beta, \phi, \eta, \alpha, \rho$

**Output**: 6-panel figure showing equilibria variation with $\eta$ and $\rho$

**How to Run**:
```bash
jupyter notebook equilibriaSHARimport.ipynb
```
Execute cells sequentially. Modify `params` variable to explore different parameter combinations.

---

### 2. `two_strain_meanfield.py`

**Purpose**: Deterministic mean-field ODE model for two-strain SHAR dynamics

**Model Structure**:
- 6 compartments: S (Susceptible), Hw (Hospitalized wild-type), Aw (Asymptomatic wild-type), Hm (Hospitalized mutant), Am (Asymptomatic mutant), R (Recovered)
- Uses `scipy.integrate.odeint` to solve the ODE system
- Integrates from $t=0$ to $t=1000$ days

**Key Parameters**:
- Recovery rates: $\gamma_w = 1/5$, $\gamma_m = 1/7$ (days⁻¹)
- Transmission rates: $\beta_w, \beta_m$ (scaled by recovery rates)
- Severity fractions: $\eta_w = 0.45$, $\eta_m = 0.2$ (fraction hospitalized)
- Mutation rate: $\epsilon = 0.0001$ (wildtype to mutant mutation probability)
- Import rate: $\rho = 10^{-5}$ (per capita per day)
- Waning immunity: $\alpha = 1/60$ (days⁻¹)

**How to Run**:
```bash
python two_strain_meanfield.py
```

**Output**: Time-series plots showing dynamics of hospitalized and asymptomatic individuals over 1000 days

**Customization**: Edit parameter values in the script to explore different scenarios

---

### 3. `two_strainG.py`

**Purpose**: Stochastic Gillespie simulation of the two-strain SHAR model

**Model Components**:
- Uses GillesPy2 `Model` class for stochastic simulation
- Implements 21 reactions modeling spillover, infection, and recovery events
- Parameters p1-p16 encode interaction rates between compartments
- Species: S, Hw, Aw, Hm, Am, R

**Solver Options**:
- `TauLeapingSolver`: Fast approximate stochastic solver (recommended for speed)
  - Generates 500 trajectories by default
  - Good balance between accuracy and computational time
- `ODESolver`: Deterministic ODE solver for reference solution

**Key Parameters**:
- Population size: $N = 100,000$
- Time span: 0-300 days (301 time points)
- Number of stochastic trajectories: 500
- Random seed: 10 (reproducible results)

**How to Run**:
```bash
python two_strainG.py
```

**Output**: 
- 500 stochastic trajectories
- Deterministic ODE solution overlay (dashed line)
- 3-panel comparison plots:
  - Hospitalized (Hw + Hm)
  - Asymptomatic (Aw + Am)
  - Total infected ((Hw + Aw) + (Hm + Am))

**Customization**: 
- Modify `num_trajectories` to change number of simulations
- Adjust solver parameters for different accuracy/speed tradeoffs
- Change `random_seed` for different random realizations

---

### 4. `AmSobol.ipynb`

**Purpose**: Global sensitivity analysis of $A_m^*$ (mutant asymptomatic equilibrium) using Sobol' indices


**How to Run**:
```bash
jupyter notebook AmSobol.ipynb
```
Execute cells sequentially. The analysis requires ~22,000 model evaluations; this may take several minutes.

**Output**: 
- Console: First-order (S1) and total-order (ST) Sobol indices with 95% confidence intervals
- Plots: 
  - Bar chart comparing first-order vs total-order indices
  - Histogram of output distribution ($A_m^*/N$)

**Notes**: 
- "High convergence" warnings can be ignored; they indicate uncertainty in index estimates but solution is typically valid
- Results quantify which parameters most influence $A_m^*$ at equilibrium

---

### 5. `Contour_plots(R0_and_equilibria).ipynb`

**Purpose**: Visualize basic reproduction number ($R_0$) and equilibria as functions of key transmission parameters

**Plots Generated**:
1. $R_0$ contour plot as a 2D function of hospitalization rate ($\eta$) and asymptomatic transmission multiplier ($\phi$)
2. Four 2×2 grids showing hospitalized ($H^*$) and asymptomatic ($A^*$) equilibria:
   - Equilibria computed for $\rho = 10^{-6}$, $10^{-5}$, $10^{-4}$, $10^{-3}$
   - Each grid shows both $H^*$ and $A^*$ as functions of $\eta$ and $\phi$

**Computational Details**:
- Mesh grid: 100 points × 100 points (eta_range × phi_range)
- Total computations: ~40,000 equilibrium evaluations

**How to Run**:
```bash
jupyter notebook "Contour_plots(R0_and_equilibria).ipynb"
```
Execute cells sequentially. Provides intuitive visualization of parameter sensitivity.

**Output**: High-resolution contour plots showing parameter landscapes

---

## Model Equations Reference

All scripts implement variations of the SHAR model (Susceptible-Hospitalized-Asymptomatic-Recovered). The basic single-strain model equations are:

$$\frac{dS}{dt} = -\beta \frac{S}{N}(H + \phi A + \rho N) + \alpha R$$

$$\frac{dH}{dt} = \eta \beta \frac{S}{N}(H + \phi A + \rho N) - \gamma H$$

$$\frac{dA}{dt} = (1-\eta) \beta \frac{S}{N}(H + \phi A + \rho N) - \gamma A$$

$$\frac{dR}{dt} = \gamma(H + A) - \alpha R$$

**Parameter Meanings**:
- $\beta$: Transmission rate from hospitalized individuals
- $\phi$: Relative transmission from asymptomatic individuals
- $\rho$: External import rate (pathogen pressure)
- $\gamma$: Recovery rate
- $\eta$: Fraction of infected becoming hospitalized
- $\alpha$: Loss of immunity rate (waning immunity)
- $N$: Total population size

The two-strain versions duplicate these equations for wild-type (w) and mutant (m) strains, including cross-strain infection terms and mutation dynamics ($\epsilon$ probability of mutation from wildtype to mutant).

---

## Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| `ModuleNotFoundError: No module named 'gillespy2'` | Run `pip install gillespy2` |
| `ModuleNotFoundError: No module named 'SALib'` | Run `pip install SALib` |
| `ModuleNotFoundError: No module named 'jupyter'` | Run `pip install jupyter` |
| Slow execution in `AmSobol.ipynb` | Normal; Sobol analysis requires ~22,000 model evaluations. Typically takes 5-10 minutes. |
| Jupyter kernel crashes with GillesPy2 | Try using `TauLeapingSolver` instead of stochastic solver for faster execution |
| `scipy.integrate.odeint` convergence warnings | Usually safe to ignore; solution is typically valid. Indicates system reached requested tolerance. |
| AttributeError with GillesPy2 solvers | Ensure GillesPy2 ≥ 3.0 is installed: `pip install --upgrade gillespy2` |

---

## Data Output and Visualization

All notebooks and scripts generate matplotlib figures. To save plots programmatically:

```python
# Save as high-quality PDF (best for publications)
plt.savefig("filename.pdf", bbox_inches='tight', dpi=300)

# Save as PNG (suitable for presentations)
plt.savefig("filename.png", dpi=300, bbox_inches='tight')

# Save as EPS (vector format compatible with LaTeX)
plt.savefig("filename.eps", bbox_inches='tight')
```

---

## Citation

If you use this code in your research, please cite:

Saldaña, F., Stollenwerk, N., & Aguiar, M. (2024). Modelling COVID-19 mutant dynamics: understanding the interplay between viral evolution and disease transmission dynamics. *Royal Society Open Science*, 11(10), 240919.

---

## License

See LICENSE file for details.





