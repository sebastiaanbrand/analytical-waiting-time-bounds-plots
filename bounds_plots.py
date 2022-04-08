import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
import probability_tools as prob_tools
from repeater_chain_analyzer import RepeaterChainCalculator

plot_format = 'pdf'
results_folder = 'results/'


def calc_or_load(pgen, pswap, n, trunc):
    outputfolder = results_folder + 'pgen{}_pswap{}_n{}/'.format(pgen,pswap,n)
    Path(outputfolder).mkdir(parents=True, exist_ok=True)
    calculator = RepeaterChainCalculator(n=n, pgen=pgen, pswap=pswap, 
                                         trunc=trunc, outputfolder=outputfolder)

    # check if data already exists
    if Path(outputfolder + 'data/pmfs.csv').exists():
        calculator.load_from_folder(outputfolder)
        # if trunc is not high enough
        if (calculator.trunc < trunc):
            calculator = RepeaterChainCalculator(n=n, pgen=pgen, pswap=pswap, 
                                         trunc=trunc, outputfolder=outputfolder)
            calculator.calculate()
    else:
        calculator.calculate()
    
    return calculator


def mu_0(pgen, discrete):
    if discrete:
        return ( 1.0 - (1.0 / (np.log(1.0 - pgen))) ) # upper bound on mu0
    else:
       return (1.0 / pgen) # lower bound on mu0

def nu_0(pgen):
    return (3.0 - 2.0*pgen) / (pgen * (2.0 - pgen))

def mean_upper_bound(pgen, pswap, n):
    return ( 3.0 / (2.0*pswap) )**n  *  mu_0(pgen, True)

def mean_lower_bound(pgen, pswap, n):
    return (1.0/pswap) * ( (3.0 - 2.0*pswap) / (pswap*(2.0-pswap)) )**(n-1) * nu_0(pgen)

def m_upper(pgen, pswap, n):
    return (3.0/2.0) * ( 3.0 / (2.0*pswap) )**(n-1)  *  mu_0(pgen, True)

def m_lower(pgen, pswap, n):
    return ( (3.0 - 2.0*pswap) / (pswap*(2.0 - pswap)) )**(n-1)  *  nu_0(pgen)

def mean_triv_upper_bound(pgen, pswap, n):
    return ( 2.0 / pswap)**n  *  (1.0 / pgen)


def mean_triv_lower_bound(pgen, pswap, n):
    return  (1.0 / pswap)**n  *  (1.0 / pgen)


def mean_numerical(pgen, pswaps, n, trunc):
    res_upper = []
    res_lower = []
    for pswap in pswaps:
        calculator = calc_or_load(pgen, pswap, n, trunc)
        res_upper.append(calculator.mean_upper_bound(n))
        res_lower.append(calculator.mean_lower_bound(n))
    
    return np.array(res_lower), np.array(res_upper)


def plot_mean_bounds():
    pgen = 0.5
    n = 4
    trunc = int(1e5)
    pswaps = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    

    scaling = 4.25
    plt.figure(figsize=(scaling, scaling*0.75))
    ub_color = 'tab:blue'
    lb_color = 'tab:orange'
    num_color = 'tab:green'
    legend = []


    # (A) plot absolute bounds
    # 1. trivial upper bound
    triv_upper_bound = mean_triv_upper_bound(pgen, pswaps, n)
    plt.plot(pswaps, triv_upper_bound, linestyle=':',  color=ub_color)
    legend.append('$1/p_{gen}(2/p_{swap})^n$')

    # 2. trivial lower bound
    triv_lower_bound = mean_triv_lower_bound(pgen, pswaps, n)
    plt.plot(pswaps, triv_lower_bound, linestyle=':',  color=lb_color)
    legend.append('$1/p_{gen}(1/p_{swap})^n$')

    # 3. new upper bound
    upper_bound      = mean_upper_bound(pgen, pswaps, n)
    plt.plot(pswaps, upper_bound,      linestyle='-.', color=ub_color)
    legend.append('new upper bound (prop ..)')

    # 4. new lower bound
    lower_bound      = mean_lower_bound(pgen, pswaps, n)
    plt.plot(pswaps, lower_bound,      linestyle='-.', color=lb_color)
    legend.append('new lower bound (prop ..)')

    # 5. numberical mean
    num_lb, num_ub   = mean_numerical(pgen, pswaps, n, trunc)
    plt.plot(pswaps, num_ub,           linestyle='-',  color=num_color)
    plt.plot(pswaps, num_lb,           linestyle='-',  color=num_color)
    legend.append('numerical mean')

    plt.ylabel('E[T]')
    plt.xlabel('$p_{swap}$')
    plt.title('N = {}, pgen = {}'.format(2**n, pgen))
    plt.legend(legend)
    plt.tight_layout()
    plt.savefig('{}mean_bounds_plot_N{}.{}'.format(results_folder,2**n, plot_format))
    plt.clf()


    # check if mean lower bound exceeds numerical upper bound
    issue_values = pswaps[np.where((num_ub - lower_bound) < 0)]
    if (len(issue_values) > 0):
        print("mean lower check fails for pswap = " + str(issue_values))

    # check if mean upper bound is lower than the numerical lower bound
    issue_values = pswaps[np.where((upper_bound - num_lb) < 0)]
    if (len(issue_values) > 0):
        print("mean upper check fails for pswap = " + str(issue_values))


    # (B) plot ratios
    legend = []
    scaling = 4.25
    plt.figure(figsize=(scaling, scaling*0.95))

    # 0. reference
    num_ratio = num_ub / num_lb
    if (np.max(num_ratio) > 1.01):
        print("WARNING: max (num upper / num lower) = {}".format(np.max(num_ratio)))
    plt.plot(pswaps, np.ones(len(pswaps)), color=num_color)
    legend.append('numerical')

    # 1. trivial upper bound
    triv_ub_ratio = triv_upper_bound / num_ub
    plt.plot(pswaps, triv_ub_ratio, linestyle=':', color=ub_color)
    legend.append('previous upper bound')

    # 2. trivial lower bound
    triv_lb_ratio = triv_lower_bound / num_ub
    plt.plot(pswaps, triv_lb_ratio, linestyle=':', color=lb_color)
    legend.append('previous lower bound')

    # 3. new upper bound
    prop_ub_ratio = upper_bound / num_ub
    plt.plot(pswaps, prop_ub_ratio, linestyle='-.', color=ub_color)
    legend.append('new upper bound')

    # 4. new lower bound
    prop_lb_ratio = lower_bound / num_ub
    plt.plot(pswaps, prop_lb_ratio, linestyle='-.', color=lb_color)
    legend.append('new lower bound')

    plt.yscale('log', base=2)
    plt.yticks(ticks=[8, 4, 2, 1, 0.5, 0.25, 0.125], 
               labels=['8', '4', '2', '1', '1/2', '1/4', '1/8'])
    plt.ylabel('ratio to numerical mean')
    plt.xlabel('$p_{swap}$')
    plt.legend(legend, loc='upper center', bbox_to_anchor=(0.5,-0.2))
    plt.subplots_adjust(bottom=0.4)
    plt.tight_layout()
    plt.savefig('{}mean_bounds_ratios_plot.{}'.format(results_folder, plot_format))
    plt.title('N = {}, pgen = {}'.format(2**n, pgen))
    plt.savefig('{}mean_bounds_ratios_plot_title.{}'.format(results_folder, plot_format))
    plt.clf()
    


def co_cdf_markov_bound(pgen, pswap, n, t, trivial=True):
    if (trivial):
        expectation = mean_triv_upper_bound(pgen, pswap, n)
    else:
        expectation = mean_upper_bound(pgen, pswap, n)
    markov_ub = (expectation / (t + 1))
    markov_ub = np.clip(markov_ub, 0, 1)
    return markov_ub

def co_cdf_upper_bound(pgen, pswap, n, t):
    m_up = m_upper(pgen, pswap, n)
    tail_ub = np.exp( pswap - ((pswap*t)/m_up) )
    tail_ub = np.clip(tail_ub, 0, 1)
    return tail_ub

def co_cdf_lower_bound(pgen, pswap, n, t):
    m_low = m_lower(pgen, pswap, n)
    tail_lb = np.exp( ((-pswap*t)/m_low) * (1./(1.-pswap)) )
    tail_lb = np.clip(tail_lb, 0, 1)
    return tail_lb


def co_cdf_numerical(pgen, pswap, n, trunc):
    calculator = calc_or_load(pgen, pswap, n, trunc)
    pmf = calculator.pmf[n]
    co_cdf = 1 - prob_tools.pmf_to_cdf(pmf)
    return co_cdf


def plot_tail_bounds():
    pgen = 0.1
    n = 4
    pswap = 0.2

    for pswap in [0.2, 0.5]:
        if pswap == 0.5:
            trunc = 3000
        elif pswap == 0.2:
            trunc = 140000  # for pswap = 0.2 (pgen=0.1, n=4), set trunc = 140000
        else:
            raise NotImplementedError(
                "Not found good value for trunc for pswap other than 0.2 and 0.5," +\
                "(at pgen = 0.1 and n=4)")

        support = np.array(range(1, trunc+1))
        plot_until = 0.01 # clip where numerical support reaches 99% of pm

        # styling stuff
        ub_style    = '-.'
        ub_color    = 'tab:blue'
        mkv_style   = ':'
        lb_style    = '-.'
        lb_color    = 'tab:orange'
        num_style   = '-'
        num_color   = 'tab:green'
        scaling = 4.25
        fig, ax = plt.subplots(1,1, figsize=(scaling, scaling*0.75))
        #fig.figure(figsize=(scaling, scaling*0.75))

        # numerial co-CDF
        co_cdf_num = co_cdf_numerical(pgen, pswap, n, trunc)[1:] # remove t = 0
        clip_at = np.sum(co_cdf_num >= plot_until)
        if (clip_at == len(support)):
            print("NOTE: only captured {:.3f} of probability mass".format(1-co_cdf_num[-1]))
        ax.plot(support[:clip_at], co_cdf_num[:clip_at], linestyle=num_style, color=num_color, label='numerical')

        # Markov's inequality (using "trival" mean ub)
        markov_ub_triv = co_cdf_markov_bound(pgen, pswap, n, support, True)
        ax.plot(support[:clip_at], markov_ub_triv[:clip_at], linestyle=(0,(1,5)), color=ub_color, label='previous upper bound')

        # Markov's inequality upper bound (using new mean ub)
        markov_ub_imp = co_cdf_markov_bound(pgen, pswap, n, support, False)
        ax.plot(support[:clip_at], markov_ub_imp[:clip_at], linestyle=mkv_style, color=ub_color, label='previous upper bound (improved)')

        # co-CDF upper bound
        co_cdf_ub = co_cdf_upper_bound(pgen, pswap, n, support)
        ax.plot(support[:clip_at], co_cdf_ub[:clip_at], linestyle=ub_style, color=ub_color, label='new upper bound')

        # co-CDF lower bound
        co_cdf_lb = co_cdf_lower_bound(pgen, pswap, n, support)
        ax.plot(support[:clip_at], co_cdf_lb[:clip_at], linestyle=lb_style, color=lb_color, label='new lower bound')

        # axis info and stuff
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,4), useMathText=True)
        ax.set_xlabel('$t$ $(L_0/c)$')
        ax.set_ylabel('$\\Pr(T > t)$')
        fig.tight_layout()
        fig.savefig('{}tail_bounds_plot_pswap{}.{}'.format(results_folder,pswap,plot_format))
        ax.set_title('N = {}, pgen = {}, pswap = {}'.format(2**n, pgen, pswap))
        fig.tight_layout()
        fig.savefig('{}tail_bounds_plot_pswap{}_title.{}'.format(results_folder,pswap,plot_format))   

        # plot legend in separate figure
        handles, labels = ax.get_legend_handles_labels()
        fig_legend = plt.figure(figsize=(3, 1.3))
        axi = fig_legend.add_subplot(111)
        fig_legend.legend(handles, labels, loc='center', scatterpoints=1)
        axi.axis('off')
        fig_legend.canvas.draw()
        fig_legend.savefig('{}tail_bounds_plot_legend.{}'.format(results_folder,plot_format))
        
        fig_legend.clf()
        fig.clf()



if __name__ == '__main__':
    # parse args
    if (len(sys.argv) <= 1):
        print("Please specify which plots to generate [mean|tail]")
        exit(1)
    else:
        which_plots = sys.argv[1]

    Path(results_folder).mkdir(parents=True, exist_ok=True)
    if (which_plots == 'mean'):
        plot_mean_bounds()
    elif (which_plots == 'tail'):
        plot_tail_bounds()
    else:
        print("The first argument should be [mean|tail]")
        exit(1)  
