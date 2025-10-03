import yfinance as yf
import numpy as np
from scipy.integrate import quad
import pandas as pd
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import fsolve

#initialisation des paramètres du modèle de Heston
kappa = 2.0    
theta = 0.04  
sigma = 0.3    
rho = -0.7     
V0 = 0.04 
r=0.0413     

    
# Extraction des données d'options pour SPX500
spx = yf.Ticker('^SPX')
expiration_date = '2025-10-31'
opt_chain = spx.option_chain(expiration_date)
spot_price = spx.history(period="1d")['Close'][0]
tau = (np.datetime64(expiration_date) - np.datetime64('today')).astype('timedelta64[D]').astype(int) / 365.0

#Nettoyage des données
calls_filtered = opt_chain.calls[(opt_chain.calls['strike'] >= spot_price - 2000) & (opt_chain.calls['strike'] <= spot_price + 2000)]
calls_filtered.to_csv('test.csv', index=False)
calls = pd.read_csv('test.csv')
calls['spread'] = calls['ask'] - calls['bid']
calls['weight'] = 1 / (calls['spread']**2)
market_prices = (calls['bid'].values + calls['ask'].values) / 2

# Fonction de densité pour le modèle de Heston
def densité_heston(phi, S, K, tau, r, kappa, theta, sigma, rho, i, V0):
    a = kappa * theta
    x=np.log(S)
    if i==1:
        u=0.5
        b=kappa  - rho * sigma
    else:
        u=-0.5
        b=kappa
    d = np.sqrt((rho * sigma * 1j * phi - b)**2 - sigma**2 * (2 * u * 1j * phi - phi**2))
    g= (b - rho * sigma * 1j * phi + d) / (b - rho * sigma * 1j * phi - d)
    tc=(1 - g * np.exp(d * tau))/(1 - g) 
    C = r * 1j * phi * tau + a / sigma**2 * ((b - rho * sigma * 1j * phi + d) * tau - 2 * np.log(tc))
    td=(1 - np.exp(d * tau)) / (1 - g * np.exp(d * tau))
    D = (b - rho * sigma * 1j * phi + d)*td/ sigma**2 
    f = np.exp(C + D * V0 + 1j * phi * x)
    y=np.exp(-1j * phi * np.log(K)) * f / (1j * phi)
    return np.real(y)

#Fonction de prix des calls dans le modèle de Heston
def prix_call_heston(S, K, tau, r, kappa, theta, sigma, rho, V0):
    integral1, _ = quad(lambda phi: densité_heston(phi, S, K, tau, r, kappa, theta, sigma, rho, 1, V0), 0, 200)
    integral2, _ = quad(lambda phi: densité_heston(phi, S, K, tau, r, kappa, theta, sigma, rho, 2, V0), 0, 200)
    P1 = 0.5 + (1/np.pi) * integral1
    P2 = 0.5 + (1/np.pi) * integral2
    call_price = S * P1 - K * np.exp(-r * tau) * P2
    return call_price

#fonction de calcul des prix des calls pour notre ensemble de test
def calcul_prix_calls_heston( kappa, theta, sigma, rho, V0):
    prix_calcules = []
    for index, row in calls.iterrows():
        prix = prix_call_heston(spot_price, row['strike'], tau, r, kappa, theta, sigma, rho, V0)
        prix_calcules.append(prix)
    return np.array(prix_calcules)

#fonction de calibration des paramètres du modèle de Heston
def residuals(params, calls, spot_price, r):
    heston_prices = calcul_prix_calls_heston(*params)
    weights = calls['weight'].values
    return weights * (heston_prices - market_prices)

init_params = [kappa, theta, sigma, rho, V0]
res = least_squares(residuals, init_params, args=(calls, spot_price, r), method='lm')
print(res.x, res.success, res.)

# Fonction pour calculer la volatilité implicite à partir des prix du modèle de Heston
#def implied_vol_call(market_prices, kappa, theta, sigma, rho, V0):
    #implied_vols = np.array([
        #fsolve(
            #lambda sigma_: calcul_prix_calls_heston(kappa, theta, sigma_, rho, V0)[i] - market_price,0.3
      #  ) for i, market_price in enumerate(market_prices)
   # ])
   # return implied_vols





# smile pour le marché et heston (calibré et non calibré)
plt.figure(figsize=(12, 8))
plt.plot(calls['strike'], calls['impliedVolatility'],'o', label='Market Implied Volatility')
#plt.plot(calls['strike'], implied_vol_call(market_prices, kappa, theta, sigma, rho, V0),'-', label='Heston Implied Volatility (not calibrated)')
#plt.plot(calls['strike'], implied_vol_call(market_prices, *res.x),'x', label='Heston Implied Volatility (calibrated)')
plt.legend()
plt.xlabel("Strike")
plt.ylabel("Implied Volatility")
plt.title("Volatility Smile")
plt.show()