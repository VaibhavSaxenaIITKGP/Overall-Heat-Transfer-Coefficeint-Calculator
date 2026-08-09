# Overall-Heat-Transfer-Coefficeint-Calculator
Overall heat transfer coefficient (U) for tubes embedded in a sand packed-bed TES. Combines Gnielinski/Dittus–Boelter (single-phase), Gungor–Winterton + Cooper (flow boiling), radial wall conduction (Incropera), and stagnant-bed sand-side resistance (Nu = 2).
# Overall Heat Transfer Coefficient for Embedded Tubes in Sand Packed-Bed TES

This repository contains a Python script that calculates the **overall heat transfer coefficient (U)** for tubes embedded in a sand packed-bed thermal energy storage (TES) system.

The model is based on a thermal resistance network and evaluates the contribution of:
- Inner convective heat transfer (water / steam)
- Tube wall conduction
- Sand-side heat transfer (stagnant packed bed)

---

## Physical Model

The overall heat transfer coefficient is calculated from the sum of three thermal resistances (per unit tube length):

\[
\frac{1}{U \cdot A} = R_i + R_{\text{wall}} + R_{\text{sand}}
\]

Where:

- \( R_i = \dfrac{1}{h_i \pi D_i} \) → Inner convection resistance  
- \( R_{\text{wall}} = \dfrac{\ln(D_o / D_i)}{2\pi k_{\text{cu}}} \) → Radial conduction through the tube wall  
- \( R_{\text{sand}} = \dfrac{1}{h_{\text{sand}} \pi D_o} \) → Sand-side resistance  

The sand-side heat transfer coefficient is calculated using the stagnant-bed conduction limit:

\[
Nu_{\text{sand}} = 2 \quad \Rightarrow \quad h_{\text{sand}} = \frac{2\, k_{\text{eff}}}{D_o}
\]

---

## Inner Heat Transfer Correlations

| Regime                    | Correlation                          | Reference |
|---------------------------|--------------------------------------|---------|
| Single-phase liquid/steam | Gnielinski (default) / Dittus–Boelter | Gnielinski (1976), Dittus & Boelter (1930) |
| Flow boiling              | Gungor–Winterton + Cooper            | Gungor & Winterton (1986), Cooper (1984) |

Fluid properties are obtained from **CoolProp** using the IAPWS-IF97 formulation for water/steam.

---

## Geometry & Operating Conditions (Default)

| Parameter              | Value          | Unit      |
|------------------------|----------------|-----------|
| Tube inner diameter    | 70             | mm        |
| Tube outer diameter    | 78             | mm        |
| Wall material          | Copper         | –         |
| Thermal conductivity   | 400            | W/m·K     |
| Effective sand k       | 1.5            | W/m·K     |
| Mass flow rate         | 0.1            | kg/s      |
| Pressure               | 10             | bar       |

Three characteristic regimes are evaluated:

1. Sensible water heating at 60 °C  
2. Flow boiling at vapor quality \( x = 0.5 \)  
3. Superheated steam at 500 °C  

---

## Installation

```bash
pip install CoolProp numpy pandas

Usage
Simply run the script:
Bashpython overall_U_sand_TES.py
Output
The script prints:

A summary table of ( h_i ), ( h_{\text{sand}} ), overall ( U ) (outer area basis), and the percentage contribution of the sand-side resistance for each regime.
A sensitivity analysis of overall ( U ) with respect to the effective thermal conductivity of the sand bed (( k_{\text{eff}} )).
A CSV file (overall_U_results.csv) containing the tabulated results.


Key Insight
For the given geometry and stagnant packed-bed assumption, the sand-side resistance dominates the overall heat transfer. Even under strong flow boiling conditions on the water side, the overall ( U ) remains largely controlled by ( h_{\text{sand}} ).
Improving the effective thermal conductivity of the sand bed (or enhancing particle movement) has a much larger impact on system performance than further increasing the inner heat transfer coefficient.

References

Dittus, F.W. & Boelter, L.M.K. (1930). Heat Transfer in Automobile Radiators of the Tubular Type.
Gnielinski, V. (1976). New Equations for Heat and Mass Transfer in Turbulent Pipe and Channel Flow.
Gungor, K.E. & Winterton, R.H.S. (1986). A general correlation for flow boiling in tubes and annuli. International Journal of Heat and Mass Transfer, 29(3), 351–358.
Cooper, M.G. (1984). Heat flow rates in saturated nucleate pool boiling – a wide-ranging examination using reduced properties.
Incropera, F.P. & DeWitt, D.P. – Fundamentals of Heat and Mass Transfer (radial conduction).
CoolProp: http://www.coolprop.org

