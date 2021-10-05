from matplotlib.pyplot import *
from numpy import *
from scipy.stats import norm, chisquare, poisson, chi2
from scipy.optimize import curve_fit


class analyse:
    def __init__(self,data,dbm = False):
        self.data = data
        if (dbm):
            self.data  = 10**(self.data/10)
        self.norm_fittting()
        
    def norm_fittting(self):


        n, bins    = histogram(self.data,100,density=True)
        centers    = 0.5*(bins[1:]+bins[:-1])
        P0         = [mean(self.data), sqrt(mean((self.data-mean(self.data))**2))]
        try:
            pars, cov           = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, P0)
            self.mu, self.sigma = pars
            self.mu_error       = sqrt(cov[0,0])
            self.sigma_error    = sqrt(cov[1,1])
        except:
            self.mu, self.sigma = P0
            self.mu_error       = 0
            self.sigma_error    = 0

        # (self.mu, self.sigma) = norm.fit(self.data)

    def chisqure(self,x,y):
        return mean((x-y)**2/y)

    def histogram(self,cut = None):
        n, bins, patches  = hist(self.data,100,density=True)
        xx_middle = 0.5 * (bins[1:] + bins[:-1])
        if (cut):
            this_cut = self.mu + cut
        y2 = norm.pdf(xx_middle, self.mu, self.sigma)

        l = plot(xx_middle, y2, 'r--', linewidth=2)
        title(r'$\mathrm{Histogram\ of\ data:}\ \mu=%.4f\pm%.4f,\ \sigma=%.4f\pm%.4f$' %(
            self.mu,
            self.mu_error,
            self.sigma,
            self.sigma_error))
        print(f"chisqure = {chisquare(n,y2)}")
        if (cut):
            plot([this_cut,this_cut], [min(n),max(n)],label = f"{cut}"+r"$\sigma$")
        tight_layout()
        show()
        return n, xx_middle, patches


if __name__ == "__main__":
    data = np.random.normal(1,1,2000)
    an = analyse(data)
    print(an.mu, an.sigma)
    print(mean(data),sqrt(mean((data-mean(data))**2)))
    an.histogram(1.645)