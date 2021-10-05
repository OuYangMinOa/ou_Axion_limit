# ou_Axion_limit
ou_Axion_limit is a package for calculate the Axion limit and analyse normal distribution

# Download
`
pip install ou-Axion-limit
`
# Usage

## Glimit

Input you experiment parameter as this way
```
from ou_Axion_limit import Glimit
g = Glimit()
g.B    = 8
g.Q    = 30000
g.T    = 2
g.total_time = 14 * 24 * 3600
g.SNR  = 5
g.f    = 4.74e9
g.beta = 2
g.cooling = 5 * 60
g.delta_w = 5000
g.delta_v = 1000
```
### g.information()
It gives you the limit and few parameters
```
g.information()

| ================== Parameter ==============
| f =  4.740e+09 (Frequency [Hz])
| B =      8.000 (Magnetic[T])
| V =  2.341e-04 (Cavity Volume [L])
| Q =  30000.000 (Q factor)
| T =      2.000 (Noise temp [k])
| t = 1209300.000 (Integration time [s])
| ================== OUR (SNR = 5) ==========
| g_a_gamma = 2.1790437845360718e-14 GeV^-1
| g_gamma   = 2.9111980151898598
| ================== KSVZ (SNR = 5) =========
| ksvz_g_a_gamma = 7.26049021733117e-15 GeV^-1
| ksvz_g_gamma   = -0.97
| ===========================================
```
### g.to_g_gamma(x)
convert G_a_gamma_gamma to G_gamma
```
print(g.to_g_gamma(1.3e-13))
7.367
```
### g.find_limit(10)
Find the integration time with given numbers of KSVZ limit (10 is means 10 times G_KSVZ ), and it will plot a figure to you
```
g.find_limit(10)

| ============== Parameter ==================
| f =  4.740e+09 (Frequency [Hz])
| B =      8.000 (Magnetic[T])
| V =  2.341e-04 (Cavity Volume [L])
| Q =  30000.000 (Q factor)
| T =      2.000 (Noise temp [k])
| ===========================================


[*] Total time :  1209600
[*] Seaching limit : 10 times KSVZ limit
[*] Founded Answer :
	 Integration :   20.21 [minutes] to reach 10 times KSVZ limit
	 Cold down   :  300.00 [seconds]
	 Total time  :  336.00 [hours]
	 Step size   :   79.00 [kHz]
	 Move rod    :  799.53 [steps] (should be a convert to a integer)
	 Total scan  :   63.16 [MHz]
```

## analyse

import some useful package

```
from matplotlib.pyplot import *
from numpy import *
from ou_Axion_limit import analyse
```

analyse a normal distribution

```
data = np.random.normal(1,1,2000)
an = analyse(data)
print(an.mu, an.sigma)
figure()
subplot(211)
title("data")
plot(data)
subplot(212)
an.histogram(1.645)
show()
```