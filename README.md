Disclosure: In December 2019, the Center for Economic and Policy Research (CEPR) contracted with the authors
to see if the numerical and statistical results of CEPR's November 2019 study could be independently verified. Any
analysis and interpretation of findings in this report express the sole views of the authors

# 2019-Bolivian-Elections
Respository for the 2019 Bolivian Elections code and data that was used in https://www.washingtonpost.com/politics/2020/02/26/bolivia-dismissed-its-october-elections-fraudulent-our-research-found-no-reason-suspect-fraud/

original report:https://jackrw.mit.edu/sites/default/files/documents/Bolivia_report-short.pdf

Included is cleaned trep and Computo, additional analyses will be added here as well.

Below is an updated plot, with two important changes from our prior work. First, we impose a strict discount rate on ballots left to be counted, only imputing on the verified ballots at the TREP acta. This works to produce a more conservative estimate of the final margin, biased against a MAS victory. Second, we plot out the effect of the Bayeisan imperfect prior on the outcomes. The prior on the plot is such that it is the percentage size of the imperfect prior as relative to the geographic component of the simulations. The reason we make use of an imperfect prior is due to the lack of covariates within Bolivia that we can make use of for a paramtric analysis. Therefore, though non-parametric sampling and simulations, we can demonstrate first, the baseline MAS-CC margin range where the prior is set to zero, and second, how a percentage unit change in the prior strength relative to the geographic component changes results. 

![simulated bayes](https://raw.githubusercontent.com/jw2249a/2019-Bolivian-Elections/master/simulation_bayes_plot.png)

The new plot demonstrates that as sampled and simulated, excluding the prior and drawing only from the information that we have on geography, the predicted margin would be about 9.9 percentage points, as calculated from the data before the halt in the TREP count. Note that to for MAS to have surpassed the ten percentage point threshold, a prior that is two percent of the size of the geographic component is all that is necessary, given our code. 

![simulated bayes01](https://raw.githubusercontent.com/jw2249a/2019-Bolivian-Elections/master/bayes_plot0.02.jpg)

Therefore, there is no denying that the election was close. However, the election outcome in favor of MAS given the OAS design would necessitate strong evidence to throw out 100 percent of the ballots due to all of them being fraudulant so as to demonstrate that fraud decisively affected the election outcome. 

Therefore, this plot can potentially aid in future research about a relative baseline as to what types of beliefs and data are necessary so as to disprove or demonstrate the impact of election fraud. 

