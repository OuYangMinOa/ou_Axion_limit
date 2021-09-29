# ou_Axion_limit
ou_Axion_limit is a package for calculate the Axion limit

# Download
`
pip install ou-Axion-limit
`
# Usage
```
from ou_Axion_limit.Glimit import Glimit
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
g.find_limit(10)
```
