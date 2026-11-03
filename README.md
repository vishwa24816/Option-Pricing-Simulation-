# Nifty Option Chain Simulator

Builds a full Nifty option chain (calls and puts across a strike range) from a spot price, time to expiry, implied volatility and a few other inputs. It does this with four independent pricing models, so you can compare them side by side:

1. Black-Scholes (closed form)
2. Binomial tree (Cox-Ross-Rubinstein)
3. Trinomial tree (Boyle)
4. Monte Carlo simulation (geometric Brownian motion)

The code is written as 10 Jupyter notebook cells that run top to bottom.

---

## Requirements

- Python 3.9+
- `numpy`, `pandas`, `scipy`, `matplotlib`
- Optional: `yfinance` (only if you want to pull live Nifty spot)

```bash
pip install numpy pandas scipy matplotlib
pip install yfinance   # optional
```

---

## Cell map

| Cell | Purpose |
|------|---------|
| 1 | Imports and display settings |
| 2 | Inputs: spot, days, IV, rates, strike grid, model settings |
| 3 | `bs_chain`: Black-Scholes pricer |
| 4 | `binomial_chain`: CRR binomial tree |
| 5 | `trinomial_chain`: Boyle trinomial tree |
| 6 | `mc_chain`: Monte Carlo pricer with antithetic variates |
| 7 | `to_chain`: turns call and put arrays into a chain table |
| 8 | Runs all four models and prints four chains |
| 9 | Compares models against Black-Scholes and checks put-call parity |
| 10 | Plots prices, straddle and model error |

Run cells in order. Cells 3 to 7 only define functions; cell 8 does the pricing.

---

## Inputs (Cell 2)

| Variable | Default | Meaning |
|----------|---------|---------|
| `SPOT` | 25000.0 | Nifty spot price |
| `DAYS` | 7 | Calendar days to expiry |
| `IV` | 0.14 | Annualised implied volatility (0.14 = 14%) |
| `R` | 0.065 | Risk-free rate, annualised |
| `Q` | 0.013 | Continuous dividend yield |
| `STEP` | 50 | Strike interval |
| `N_STRIKES` | 20 | Number of strikes on each side of ATM (41 strikes total) |
| `TREE_STEPS` | 500 | Time steps for both trees |
| `MC_PATHS` | 200000 | Monte Carlo paths (must be even, since antithetic pairs are used) |
| `SEED` | 42 | Random seed for reproducible Monte Carlo |

Derived values:

- `T = DAYS / 365`
- `ATM = round(SPOT / STEP) * STEP`
- `STRIKES = ATM + STEP * [-N_STRIKES ... +N_STRIKES]`

The defaults are placeholders. Set them to the current Nifty level, the expiry you care about, and the IV you see on your broker's chain.

---

## Models

All four models price European options and use the same cost-of-carry drift, `r - q`. Since Nifty index options are European-style, this matches the contract.

### 1. Black-Scholes

```
d1 = [ln(S/K) + (r - q + sigma^2/2) T] / (sigma sqrt(T))
d2 = d1 - sigma sqrt(T)

Call = S e^(-qT) N(d1) - K e^(-rT) N(d2)
Put  = K e^(-rT) N(-d2) - S e^(-qT) N(-d1)
```

Exact for the assumptions (constant vol and rates, lognormal prices). It is the benchmark the other three are compared to.

### 2. Binomial tree (CRR)

```
dt = T / N
u  = e^(sigma sqrt(dt)),  d = 1/u
p  = (e^((r-q) dt) - d) / (u - d)
```

Terminal prices are `S * u^(N-j) * d^j` for `j = 0..N`. Option values are computed at expiry and discounted backwards one step at a time:

```
V = e^(-r dt) * (p * V_up + (1-p) * V_down)
```

Recombining, so there are N+1 terminal nodes.

### 3. Trinomial tree (Boyle)

Works in log-price space with three branches (up, middle, down):

```
dx = sigma sqrt(3 dt)
nu = r - q - sigma^2/2
a  = (sigma^2 dt + nu^2 dt^2) / dx^2
b  = nu dt / dx

pu = (a + b) / 2
pd = (a - b) / 2
pm = 1 - pu - pd
```

There are 2N+1 terminal nodes at prices `S * e^(j dx)` for `j = -N..N`. Backward induction:

```
V = e^(-r dt) * (pu * V_up + pm * V_mid + pd * V_down)
```

The extra middle branch usually gives smoother convergence than the binomial tree for the same N.

### 4. Monte Carlo

Terminal prices are sampled directly from the GBM solution:

```
S_T = S exp[(r - q - sigma^2/2) T + sigma sqrt(T) Z],   Z ~ N(0,1)
```

Price is the discounted mean payoff.

- **Antithetic variates:** for every `Z` the code also uses `-Z`, which cuts variance at no extra random-number cost.
- **Common random numbers:** the same paths price every strike, so the chain is internally consistent and has no strike-to-strike noise from different samples.
- **Standard error:** `CE_stderr` and `PE_stderr` columns show the sampling error of each Monte Carlo price. Treat the standard error as approximate: it is computed over the individual paths, not over antithetic pairs.

---

## Output

### Cell 8: four chain tables

Each table is indexed by strike:

| Column | Meaning |
|--------|---------|
| `CE` | Model call price |
| `PE` | Model put price |
| `Straddle` | `CE + PE` |
| `CE_stderr`, `PE_stderr` | Monte Carlo only |

### Cell 9: validation

- **Error vs Black-Scholes:** per-strike differences for the binomial, trinomial and Monte Carlo prices, plus their maximum absolute differences.
- **Put-call parity residual:** checks `C - P = S e^(-qT) - K e^(-rT)` for each model. Black-Scholes and both trees should be near zero. Monte Carlo will show a small residual from sampling noise.

### Cell 10: plots

Four panels: call prices, put prices, straddle prices, and call error vs Black-Scholes. A dashed vertical line marks spot.

### Sample results

With the defaults (spot 25000, 7 days, IV 14%, 500 tree steps, 200k paths):

| Model | Max abs difference vs Black-Scholes |
|-------|------------------------------------|
| Binomial | about 0.10 |
| Trinomial | about 0.10 |
| Monte Carlo | about 0.6 |

Your numbers will differ if you change the inputs or the seed.

---

## Tuning accuracy and speed

- **Trees:** raise `TREE_STEPS` for tighter agreement with Black-Scholes. Runtime grows roughly linearly with steps. Trinomial memory grows with `(2N+1) x strikes`.
- **Monte Carlo:** error shrinks with the square root of paths, so 4x the paths halves the error. The simulation holds a `paths x strikes` array in memory (several are alive at once), so if memory is tight reduce `MC_PATHS` or `N_STRIKES`.
- **Reproducibility:** keep `SEED` fixed to get identical Monte Carlo output between runs.

---

## Assumptions and limitations

- **Flat volatility:** one IV is used for all strikes, so there is no skew or smile. Real Nifty chains have skew, so far OTM puts will be priced lower here than on the exchange.
- **European exercise only:** no early exercise is modelled.
- **Constant inputs:** rate, dividend yield and vol are constant over the life of the option. Dividends are a continuous yield, not discrete payments.
- **Calendar-day time:** `T = DAYS / 365`. If you prefer trading-day time, set `T` yourself in cell 2.
- **No market frictions:** no bid-ask spread, taxes, brokerage or margin.
- **Model prices, not quotes:** output will not match live LTPs exactly, since the market embeds skew, supply and demand, and expiry-day effects.
- **Very short expiries:** with tiny `T` (for example expiry day), the tree probabilities and Monte Carlo variance can behave poorly. Sanity-check against Black-Scholes.

---

## Extending it

- **Volatility smile:** pass a strike-dependent sigma. Black-Scholes and Monte Carlo are easy to adapt. The trees are vectorised across strikes with a single sigma, so a smile means running the tree once per strike.
- **American options:** in each tree's backward loop, take `max(continuation value, intrinsic value)` at every node.
- **Greeks:** use finite differences (bump spot, vol or time and reprice) for the trees and Monte Carlo, or analytic formulas for Black-Scholes.
- **Implied vol:** invert Black-Scholes numerically (for example with `scipy.optimize.brentq`) to back out IV from a market price.
- **Strategy payoffs:** the chain DataFrames can feed multi-leg payoff or P&L calculations.

---

## Troubleshooting

| Problem | Likely cause and fix |
|---------|----------------------|
| `NameError: STRIKES` | Cell 2 was not run before cell 7 or 8. |
| `display` not defined | Running outside Jupyter. Use `print(df)` instead. |
| Monte Carlo is slow or uses lots of RAM | Reduce `MC_PATHS` or `N_STRIKES`. |
| Trees disagree with Black-Scholes | Increase `TREE_STEPS`. |
| Live spot fetch fails | `yfinance` needs internet access and may be rate-limited. Enter `SPOT` manually. |

---

## Disclaimer

This is an educational pricing tool. Outputs are theoretical model values, not trading advice or live market prices.
