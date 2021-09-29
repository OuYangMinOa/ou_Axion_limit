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
        self.calculate()

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


    def find_limit(self,target):
        
        print("|","="*18,"Parameter","="*14)
        print(f"| f = {self.f:.5e} (Frequency [Hz])")
        print(f"| B = {self.B:.5f} (Magnetic[T])")
        print(f"| V = {self.V:.5e} (Cavity Volume [L])")
        print(f"| Q = {self.Q:.5f} (Q factor)")
        print(f"| T = {self.T:.5f} (Noise temp [k])")
        print("\n")
        print("[*] Total time : ",self.total_time)
        print("[*] Seaching limit :",target ,"times KSVZ limit")

        def start_exp(scan_shfit):
            spectrum    = zeros(int((scan_range + scan_windos)//grid_unit))
            weight      = zeros(int((scan_range + scan_windos)//grid_unit))
            sigma       = zeros(int((scan_range + scan_windos)//grid_unit))
            shift_index = int(scan_shfit/grid_unit)
            start_index = 0
            
            f_now = f_start
            noise        = random.normal(0,1,windo_index)
            while f_now <= f_end:
                noise        = np.zeros(windo_index) +1
                scan_f       = linspace(f_now - scan_windos/2, f_now + scan_windos/2, windo_index)
                sigma_this   = 1
                
                w   = Lorentz(scan_f, f_now, self.Q) / sigma_this**2
                noise *= w
                sig = sigma_this**2 * (w)**2
                
                spectrum[start_index:start_index+windo_index] = spectrum[start_index:start_index+windo_index] + noise*w
                weight[  start_index:start_index+windo_index]   = weight[start_index:start_index+windo_index] + w
                sigma[   start_index:start_index+windo_index]    = sigma[start_index:start_index+windo_index] + sig
                
                f_now       += scan_shfit
                start_index += shift_index
            return linspace(f_start, f_end, int((scan_range + scan_windos)//grid_unit)), spectrum, sigma

        def Lorentz(x,fr,Q):
            return 1 / (1 + ( 2*Q*(x/fr-1))**2)
        def merging(x,y,merge):
            out_x, out_y = [], []
            for i in range(len(x)//merge):
                out_x.append(mean(x[i*merge:(i+1)*merge]))
                out_y.append(sum( y[i*merge:(i+1)*merge]))
            return array(out_x), array(out_y)

        grid_unit   = 1e3
        f_start     = self.f - 200e3*100
        f_end       = self.f + 200e3*100
        scan_windos = 2e6
        windo_index = int(scan_windos/grid_unit)
        scan_range  = f_end - f_start
        start_index = int((scan_windos/2) / grid_unit)
        end_index   = start_index + int(scan_range / grid_unit)

        x,spec, sigm = start_exp(self.f / self.Q / 2)
        # print(spec,sigm)
        
        _, spec = merging(x,spec,5)
        x, sigm = merging(x,sigm,5)
        y = np.divide(spec, sqrt(sigm), out=np.zeros_like(spec), where=sqrt(sigm)!=0)

        self.t = 600*60
        self.cal_single_point()
        with np.errstate(divide='ignore'):
            g_a_gamma = sqrt(self.sigma*5/self.shift/y) * 1e9

        g_gamma = (pi * self.big_A * self.big_A) * g_a_gamma / self.ma / self.alpha   * 1e-9
        # print(min(g_gamma/0.97))

        new_t = self.t * ( min(g_gamma/0.97) / target)  **4
        print(f"[*] Founded Answer :")
        print(f"\t Integration : {(new_t/60):7.2f} [minutes] to reach {target} times KSVZ limit")
        print(f"\t Cold down   : {self.cooling:7.2f} [seconds]")
        print(f"\t Total time  : {self.total_time/3600:7.2f} [hours]")
        print(f"\t Step size   : {self.f / self.Q / 2*1e-3:7.2f} [kHz]")
        print(f"\t Move rod    : {self.total_time / (new_t + self.cooling):7.2f} [steps]")
        print(f"\t Total scan  : {self.total_time / (new_t + self.cooling) * self.f / self.Q / 2*1e-6:7.2f} [MHz]")


        self.t = new_t

        assert new_t < self.total_time ,"Integration time is larger then total time"


        f_start     = self.f - self.f / self.Q / 2*(self.total_time / (new_t + self.cooling)//2)
        f_end       = self.f + self.f / self.Q / 2*(self.total_time / (new_t + self.cooling)//2)
        scan_range  = f_end - f_start
        start_index = int((scan_windos/2) / grid_unit)
        end_index   = start_index + int(scan_range / grid_unit)
        x,spec, sigm = start_exp(self.f / self.Q / 2)
        _, spec = merging(x,spec,5)
        x, sigm = merging(x,sigm,5)
        y = np.divide(spec, sqrt(sigm), out=np.zeros_like(spec), where=sqrt(sigm)!=0)

        self.cal_single_point()
        with np.errstate(divide='ignore'):
            g_a_gamma = sqrt(self.sigma*5/self.shift/y) * 1e9
        g_gamma = (pi * self.big_A * self.big_A) * g_a_gamma / self.ma / self.alpha   * 1e-9

        figure()
        max_ratio = min(g_gamma/0.97)
        title(f"Best value : {max_ratio:.3f} ")
        plot(x*1e-9,g_gamma/0.97,label="This experiment")

        plot(x*1e-9,abs(g.ksvz_g_gamma +x*0)/0.97,"b--",label="KSVZ")
        plot(x*1e-9,abs(g.dfsz_g_gamma +x*0)/0.97,"r--",label="DFSZ")
        upper = 4
        down  = abs(g.dfsz_g_gamma +x*0)/abs(g.ksvz_g_gamma) / 4
        fill_between(x*1e-9,upper,down,color="yellow",label="model region")

        xlabel(f"Freq [GHz]")
        ylabel(r"$\frac{G_\gamma}{G_{KSVZ}}$",size=20)
        xlim(min(x*1e-9),max(x*1e-9))
        tight_layout()
        grid()
        gcf().autofmt_xdate()
        legend()
        tight_layout()
        show()


if __name__ == "__main__":
    g = Glimit()
    g.B    = 8
    g.Q    = 30000
    g.T    = 2
    g.total_time = 14 * 24 * 3600
    g.SNR  = 5
    g.f    = 4.74e9
    g.beta = 2
    g.cooling = 5 * 60
    g.Scanran = 2e6
    g.Scanwin = 2e6
    g.Scanshi = 2e6
    g.delta_w = 5000
    g.delta_v = 1000
    g.find_limit(4)