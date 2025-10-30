#  EA_P1: Rough Hawkes Heston Models  

##  Project Description  
This repository contains all the code and theoretical developments produced as part of **EA P1**, dedicated to the study, analysis, and implementation of *Rough Hawkes Heston models* for financial volatility modeling.  

This project is **currently a work in progress**.

While the implementations of the Heston and Lifted Heston models are complete, the development of the Rough Hawkes Heston model — particularly the simulation of Hawkes and rough Hawkes processes — is still ongoing.

Future updates will include:

- Full Monte Carlo simulation of the Rough Hawkes Heston dynamics;

- Implementation and testing of Hawkes process generation schemes (Ogata’s thinning and branching representations);

---

##  Abstract  
This repository focuses on the implementation and numerical simulation of the **Rough Hawkes Heston model**, a stochastic volatility framework that combines *rough volatility* and *self-exciting jumps* to capture the joint dynamics of equity indices (e.g., S&P 500) and their volatility indices (e.g., VIX).  

The rise of rough volatility models—since the seminal work of **Jim Gatheral, Thibault Jaisson, and Mathieu Rosenbaum (2018)** in [*“Volatility is Rough”*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3014764)—has demonstrated that volatility behaves like a **fractional process** with Hurst parameter \( H < 0.5 \). This provides a more realistic description of volatility clustering and long-memory effects observed in financial markets.  

Building on this paradigm, the **Rough Hawkes Heston model** introduced by **Bondi, Pulido, and Scotti (2022)** in [*“The Rough Hawkes Heston Stochastic Volatility Model”*](https://arxiv.org/abs/2210.12393) extends the classical Heston framework by incorporating both *rough behaviour* and a *Hawkes-type jump intensity*, generating endogenous and clustered jumps in volatility and returns.  

In their paper, the authors derive **semi-closed Fourier–Laplace pricing formulas** for call and put options on both the underlying and the VIX.  

The purpose of this project is to extend those results by:  
- **Developing Monte Carlo simulation methods** for the Rough Hawkes Heston dynamics to price complex products such as path-dependent or volatility-linked derivatives;  
- **Adapting Hawkes process simulation techniques** (e.g., Ogata’s thinning algorithm, branching representations);  
- **Validating and reproducing** VIX option prices obtained analytically through Monte Carlo methods.  

**Keywords:** rough volatility, Hawkes processes, stochastic volatility, self-exciting jumps, option pricing, Monte Carlo simulation.  

---

##  Project Structure  

To build up intuition and gradually address model complexity, we proceed as follows:  

1. **Heston Model** — baseline affine stochastic volatility model.  
   - Analytical pricing via Fourier transform (Heston 1993).  
   - Calibration to market data.  
   - Path simulation via Euler and exact schemes.  
   - Reference: [Steven L. Heston, *“A Closed-Form Solution for Options with Stochastic Volatility”* (1993)](https://doi.org/10.1093/rfs/6.2.327).  

2. **Lifted Heston Model** — multi-factor approximation of the rough Heston kernel.  
   - Implementation of the lifted representation and comparison with fractional kernel behavior.  
   - Reference: [El Euch & Rosenbaum, *“The Characteristic Function of Rough Heston Models”* (2019)](https://arxiv.org/abs/1609.02108).  

3. **Rough Hawkes Heston Model** — main focus of the project.  
   - Simulation of the coupled volatility–price system.  
   - Monte Carlo pricing of vanilla and exotic derivatives.  
   - Comparison between analytical Fourier-based and numerical Monte Carlo results.  
   - Reference: [Bondi, Pulido & Scotti (2022)](https://arxiv.org/abs/2210.12393).  


---
## Autors 

- Christel Astride Mallo P. — Applied Mathematics, École Polytechnique ( astride.mallopoundi@gmail.com)

- Achley Tchouala M. — Applied Mathematics, École Polytechnique (achley.tchouala-momene@polytechnique.edu)

- Felix Albert Essama O. — Applied Mathematics, École Polytechnique (essama.felix-albert@polytechnique.edu)

---

##  References  

- **Gatheral, J., Jaisson, T., & Rosenbaum, M. (2018).**  
  *Volatility is Rough.* **Quantitative Finance**, 18(6), 933–949.  
  [SSRN Link](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3014764)  

- **Bondi, A., Pulido, S., & Scotti, S. (2022).**  
  *The Rough Hawkes Heston Stochastic Volatility Model.*  
  [arXiv:2210.12393](https://arxiv.org/abs/2210.12393)  

- **Heston, S. (1993).**  
  *A Closed-Form Solution for Options with Stochastic Volatility with Applications to Bond and Currency Options.*  
  **Review of Financial Studies**, 6(2), 327–343.  
  [DOI:10.1093/rfs/6.2.327](https://doi.org/10.1093/rfs/6.2.327)  

- **El Euch, O., & Rosenbaum, M. (2019).**  
  *The Characteristic Function of Rough Heston Models.*  
  **Mathematical Finance**, 29(1), 3–38.  
  [arXiv:1609.02108](https://arxiv.org/abs/1609.02108)

---

##  Architecture  

The `Option` class is the **base abstract class** from which all specific model implementations (Heston, Lifted Heston, Rough Hawkes Heston) inherit.  

```python
"""
Abstract class Option.

Parameters:
    phi : +1 for Call, -1 for Put
    t0  : evaluation date
    T   : maturity
    S0  : spot price at t0
    K   : strike price
    r   : risk-free rate
    q   : dividend yield
    market_price : observed market price at t0
    vega : observed market vega at t0
"""
