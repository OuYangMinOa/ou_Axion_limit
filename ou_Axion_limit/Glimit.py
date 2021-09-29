from matplotlib.pyplot import *
from numpy import *



class Glimit:
    def __init__(self):
        self.h      = 6.626e-34 # J/k
        self.h_bar  = 4.135e-15/(2*pi) # ev
        self.k      = 1.38e-23 # J
        self.rho    = 0.45e15   # axion density (GeV/cc)
        self.big_A  = 78e6  #   78Mev
        self.C      = 3e8    # The Speed of light
        self.alpha  = 1/137  # fine structure const
        self.mu     = 4*pi * 1e-7 
        self.e      = 1.6e-19
        self.he     = self.h/self.e
        self.cm     = 0.6527    # one hour (s)

        self.s          = 0      # S₁₁
        self.T          = .3     # physical + noise (k) 
        self.f          = 5e9    # 5Ghz Frequency (10⁹ Hz)
        self.B          = 8      # Magnetic (T)
        self.V          = 120*pi*621e-9 #48338e-9 # Volume (L)
        self.Q          = 25000 # Q factor 
        self.delta_v    = 1000   # 5kHz
        self.delta_w    = 5000   # 5kHz
        self.cooling    = 0
        self.SNR        = 5
        self.beta       = 1
        self.total_time = 86400 * 7
        self.Scanwin    = 2e6
        self.Scanran    = 2e6
        self.Scanshi    = 2e6
        
        self.ksvz_g_gamma = -0.97
        self.dfsz_g_gamma = 0.36

    def calculate(self):
        
        self.w         = 2*pi*(self.f)
        self.ma        = self.h_bar * self.w
        self.t         = self.total_time / self.Scanran * self.Scanshi - self.cooling

        assert self.t  > 0,"No time for cooling down"

        self.g_KSVZ    = 0.97 * self.ma * self.alpha /(pi * self.big_A * self.big_A)
        self.Na        = self.delta_v * self.t
        self.Ns        = self.Scanwin/self.Scanshi
        self.sigma     = self.k * self.T * self.delta_w/sqrt(self.Ns * self.Na) 
        
        self.shift =  ((self.h_bar*self.C)**3)*self.rho / (self.ma**2)  * (1/self.mu) * (self.B**2) * \
        self.V * self.w* self.cm * self.Q *self.beta / (1+self.beta)
        
    def cal_single_point(self):
        self.Na    = self.delta_v * self.t
        self.Ns    = 1
        self.sigma = self.k * self.T * self.delta_w/sqrt(self.Na) 
        self.shift =  ((self.h_bar*self.C)**3)*self.rho / (self.ma**2)  * (1/self.mu) * (self.B**2) * \
        self.V * self.w* self.cm * self.Q *self.beta / (1+self.beta)

    def g_a_gamma(self):
        return sqrt( self.SNR*self.sigma/self.shift) * 1e9
    
    def g_gamma(self):
        return  (pi * self.big_A * self.big_A) * self.g_a_gamma() / self.ma / self.alpha   * 1e-9    
    
    def ksvz_g_a_gamma(self):
        return 0.97 * self.ma * self.alpha /(pi * self.big_A * self.big_A)  * 1e9
        
    def information(self):
        self.calculate()
        print("|","="*18,"Parameter","="*14)
        print(f"| f = {self.f:.5e} (Frequency [Hz])")
        print(f"| B = {self.B:.5f} (Magnetic[T])")
        print(f"| V = {self.V:.5e} (Cavity Volume [L])")
        print(f"| Q = {self.Q:.5f} (Q factor)")
        print(f"| T = {self.T:.5f} (Noise temp [k])\n|")
        print("|","="*18,f"OUR (SNR = {self.SNR:d})","="*10)
        print(f"| g_a_gamma = {self.g_a_gamma() } GeV^-1")
        print(f"| g_gamma   = {self.g_gamma()}\n|")
        print("|","="*18,f"KSVZ (SNR = {self.SNR:d})","="*9)
        print("|",f"ksvz_g_a_gamma = {self.ksvz_g_a_gamma()} GeV^-1")
        print("|",f"ksvz_g_gamma   = {self.ksvz_g_gamma}\n|")
        print("|","="*43)
