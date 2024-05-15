## Selection on unobservables designs

#### Key models

1. Random Effect and Fixed Effect Models (RE & FE)
2. Difference-in-Differences Models (DID)
3. Synthetic Controls
4. Instrumental Variable (IV)
5. Regression Discontinuity Design (RDD)



### 1 Random Effect and Fixed Effect Models

#### Notations

$i$: units/individuals

$t$: time

$Y_{it}$: outcome of unit $i$ in period $t$

$X_{it}$: a vector of explanatory variables of unit $i$ in period $t$

- Key distinction between fixed effects and random effects：
  - Fixed effect: we model the correlation between the individual effects $c_i$ and the covariates $X_{it}$
  - Random effect: we assume the individual effects $c_i$ and the covariates $X_{it}$​ are independent

#### Random Effect

- **Key Assumption**s: 

1. Strict Exogeneity: $\mathbb{E}[\varepsilon_{it}|X_{i1},...,X_{iT},c_i]=0.$​
2. Uncorrelated Effects: $\mathbb{E}[c_i|X_{1i},...,X_{iT}]=0$.

```stata
reghdfe y treat x1 x2 x3, noabsorb vce(cluster country)
est store m1
```

#### Fixed Effect

- **Key Assumption**: 

1. Strict Exogeneity: $\mathbb{E}[\varepsilon_{it}|X_{i1},...,X_{iT},c_i]=0.$​

But the fixed effect $c_i$ and covariates $X_{1i},...,X_{iT}$ are correlated.

The idea behind fixed effects is that we want to either estimate the $c_i$​ parameters (so that we can control for them) or just get rid of them altogether (so that we don't have to worry about them).

$Y_{it}=\beta D_{it}+c_i+\tau_t+\varepsilon_{it}$

- Within Estimator: use within-individual variation and demeaning transformation to produce the same FE estimator
- Between Estimator: estimate $\beta$ using only between-individual variation. It can be implemented by running a regression on individual's averages.

```stata
reghdfe y treat x1 x2 x3 if age>=50 & age<=70, a(county year) vce(cl county)
//you can also use vce(robust)
est store m //save the model

//Or a more comprehensive note
reghdfe ln_fata primary, absorb(state) vce(robust) //state fixed effect
reghdfe ln_fata primary, absorb(year) vce(robust) //year fixed effect
reghdfe ln_fata primary, absorb(state year) vce(robust) //state + year fixed effect
```



### 2 Differences-in-Difference

**Notation:** 

$i$: individual

$c$: city

$t$: time

We write $Y_{ict}=\bar{Y}+\tau D_{ct}+\gamma_c+\delta_t+\varepsilon_{ct}+u_{ict}$​

- The inclusion of $\varepsilon_{ct}$ guarantees that $\bar{u}_{ct}=0$
- We want to see how much change we would have expected in the treated group if no treatment had occurred:  (treated group after-treated group before)-(untreated group after-untreated group before)=$(\bar{Y}_{11}-\bar{Y}_{10})-(\bar{Y}_{01}-\bar{Y}_{00})$, $\bar{Y}_{ct}=\frac{1}{N_{ct}}\sum_i Y_{ict}$​
- Compared with FE, $D_i$ is changed to an interaction term between `treated` and `after`

**Key Assumption**: the parallel trends $\mathbb{E}[\varepsilon_{11}-\varepsilon_{10}]=\mathbb{E}[\varepsilon_{01}-\varepsilon_{00}]$, 

- saying that if no treatment had occurred, the difference between the treated group and the untreated group would have stayed the same in the post-treatment period as it was in the pre-treatment period/ the treated group and untreated groups had similar trajectories for the dependent variable before treatment

- This assumption is inherently unobservable. It's about the counterfactual of what would have happened if treatment had not occurred.

- How to evaluate the parallel trends assumption? Event-study!

  【Also, even if we show that the trends are parallel between the two groups before the treatment, that does not necessarily mean that the parallel trends assumption will hold, $\because$ we do not observe the counter-factual trends post treatment in most empirical setting.】

  We should not say "the assumption is valid" but say "we are assessing the plausibility of the assumption"

  $y_{it}=\sum_{k\geq maxLeadTime, k\not=1}^{k=maxLagTime}\delta_kD_{it}^k+\mu_i+\rho_t+\varepsilon_{it}$​

  - Define $s_i$ as the year when city $i$ was first assigned the policy
  - $D_{it}^{-maxLeadTime}=1$ if $t-s_i\leq -maxLeadTime$ and $0$ otherwise
  - $D_{it}^{k}=1$ if $t-s_i=k$ and $0$ otherwise, where $k=-maxLeadTime+1, ..., maxLagTime-1$
  - $D_{it}^{maxLagTime}=1$ if $t-s_i\geq maxLagTime$ and $0$ otherwise
  - $\mu_i$ is the city-level fixed effect
  - $\rho_t$ is the year-level fixed effect

  ```stata
  //B2-Event Study-to test the parallel trend
  gen n=year if primary==1
  bysort state:egen event=min(n) //get the treatment time
  gen TimeToTreat=year-event
  //Three approahes to evaluate the parallel trend
  ////Approach 1: hdfe
  eventdd ln_fata, timevar(TimeToTreat) method(hdfe, absorb(state year) cluster(state)) accum lags(10) leads(10) graph_op(ytitle("Event study of fatality rates")) 
  ////Approach 2: fe
  eventdd ln_fata i.year, timevar(TimeToTreat) method(fe, cluster(state)) accum lags(10) leads(10) graph_op(ytitle("Event study of fatality rates")) 
  ////Approach 3: ols
  eventdd ln_fata i.year i.state, timevar(TimeToTreat) method(, cluster(state)) accum lags(10) leads(10) graph_op(ytitle("Event study of fatality rates")) 
  ```

**Idea of DID**:  use the change in the untreated group to represent all non-treatment changes in the treated group.

Implementation:

1. The classic approach to estimate DID is the two-way fixed effects

$Y_{ict}=\alpha+\tau D_{ct}+\gamma 1(c=1)+\delta 1(t=1)+\varepsilon_{ct}+u_{ict}$

- $\tau$​ is the treatment effect, which tells us how much bigger the treated group effect is in the after-treatment than in the before-period

- remember to cluster the standard error

- ```stata
  egen stateid=group(state)
  gen treated=(state==99)
  gen after=(year>=1986)
  gen did=treated*after
  
  reghdfe ln_fata did X, absorb(state year) cluster(state)
  ```

  

### 3 Synthetic Control

#### Setup

1. There is a policy that goes into effect at a particular time, but only for a particular group.
2. You use data from the pre-treatment period to adjust for differences between the treatment and control groups, and then see how they differ after treatment goes into effect.

#### **Notations**

- Suppose that we observe $J+1$ units in periods $1,2,...,T$​.

- Only unit $1$ was exposed to the treatment during periods $T_0+1,...,T$.

- The remaining $J+1$ units are untreated for all periods.

- **Our target**: to construct a synthetic control group outcome for the treated unit out of the $J$​ potential control units.

- $Y_{it}^I$: outcome that would be observed for unit $i$ at time $t$ if unit $i$ is exposed to the intervention in the periods $T_0+1$ to $T$.

- $Y_{it}^N$: potential outcome that would be observed for unit $i$ at time $t$ in the absence of the intervention.

- The treatment effect of unit $1$ at time $t$ for $t>T_0$:

  $\tau_{1t}=Y_{1t}^I-Y_{1t}^N=Y_{1t}-Y_{1t}^N$, where $Y_{1t}$ is the outcome for unit one at time $t$.

- We want to construct the counterfactual as the weighted average of the outcomes from the potential controls

  $Y_{1t}^N=\sum_{j=2}^{J+1}w_j^{*}Y_{jt}$ and then $\tau_{1t}=Y_{1t}-\sum_{j=2}^{J+1}w_j^{*}Y_{jt}$, where $W=(w_2,...,w_{J+1})^{T}$ is a set of weights, with $w_j\geq 0$ for $j=2,...,J+1$ and $\sum_{j=2}^{J+1}w_j=1$.

- useful commend: [synth2](https://www.lianxh.cn/details/1133.html) ([another useful slides](https://www.stata.com/meeting/china22-Uone-Tech/slides/China22_Guanpeng.pdf)), 

```stata
//Construct the synthetic weight
synth2 ln_fata college population unemploy beer precip snow32 totalvmt ln_fata(1982) ln_fata(1983) ln_fata(1984) ln_fata(1985), trunit(49) trperiod(1986) xperiod(1982(1)1985) nested allopt
```

- Permutation Test//Placebo Test

  ```stata
  synth2 ln_fata college population unemploy beer precip snow32 totalvmt ln_fata(1982) ln_fata(1983) ln_fata(1984) ln_fata(1985), trunit(49) trperiod(1986) xperiod(1982(1)1985) placebo(unit cut(2)) sigf(6)
  ```



### 4 Instrumental Variables Desgins

**Ideal**: to find some subset of the variation in the treatment, call it $Z_i$, that is uncorrelated with $\varepsilon_i$ (i.e., as good as randomly assigned).

**Assumptions**:

1. Exogenous variation $Cov(z_i,d_i)\not=0$: $z_i$ captures some of the variation in $d_i$. If it doesn't, then it will be of no use to us in the estimating the effect of $d_i$ on $y_i$.
2. Exclusion Restriction $Cov(z_i,\varepsilon_i)=0:$ $z_i$ is uncorrelated with $\varepsilon_i$. 

**2SLS**

- This is the most popular way to implement the IV estimator via a two-stage procedure.
- If we have 1 IV and 1 variable that we want to instrument for, 2SLS=IV, i.e., IV is a special case of 2SLS.

**Steps for the 2SLS**

- The first stage: $d_i=\gamma_1 z_i+\bold{x}_i\bold{\gamma}_2+u_i$. Then we take the predicted values of $d_i$ as $\hat{d}_i$ ($\hat{d}_i=\hat{\gamma}_1z_i+\bold{x}_i\hat{\bold{\gamma}}_2$​).
- The second stage: $y_i=\beta_0+\beta_1\hat{d}_i+\bold{x}_i\bold{\beta}_2+\varepsilon_i$, $\beta_1$​ is what we want.
- Both the first and second stages always contain the same sets of covariates.

**The reduced Form**

- this measures how $y_i$ changes as we change the instrument $z_i$, $y_i=\pi_1z_i+\bold{x}_i\bold{\pi}_2+v_i$​
- We cannot directly use $\pi_i$ to represent the treatment effect, but we can rescale it by the first stage coefficient $\hat{\beta}_{1,IV}=\frac{\hat{\pi}_1}{\hat{\gamma}_1}$.

```stata
//generate IV
gen d_hrs_285 = cond(hrs_82 > 28.5, 1, 0)
//first stage
reghdfe treatment iv x1 x2, absorb(state year) vce(robust)
est store m_iv1
//second stage
ivreghdfe outcome x1 x2 (treatment = iv), absorb(state year) vce(robust)
est store m_iv2

esttab m_iv1 m_iv2 using "iv.rtf", b(3) se(3) stat(N r2 r2_a F) nonumbers mtitles("First Stage" "Second Stage")
```



### 5 Regression Discontinuity Designs

**Notations:**

- $D_i$：a binary treatment
- $Y_i$: Outcome
- $Y_i(0)$: potential untreated outcome
- $Y_i(1)$:potential treated outcome

**Key assumption**: discontinuity assumption, the relationship between the running variable $X_i$ and the outcome $Y_i$ is smooth, i.e., $Y_i(0)$ and $Y_i(1)$ do not jump discontinuously as $X_i$​​​ changes.

- This is a compensation for the failure of the overlap assumption in the RDD.

- There are three types of RD plots to check the continuity assumption 

  1.Plot the outcomes & treatment by the running variable
  2.Plot the covariates by the running varaible
  3.Plot the density of the running variable

```stata
//-1. Outcome & Running variable-//
rdplot lnmdvalhs0_nbr hrs_82, c(28.5) p(1) graph_options(title(Outcome and Running Variable)) addlabopts(mlabsize(3))
graph save rd11.png,  replace
rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(1)

rdplot lnmdvalhs0_nbr hrs_82, c(28.5) p(2) graph_options(title(Outcome and Running Variable)) addlabopts(mlabsize(3))
graph save rd12,  replace
rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(2)

//-1. Treatment & Running variable-//
rdplot npl2000 hrs_82, c(28.5) p(1) graph_options(title(Treatment and Running Variable)) addlabopts(mlabsize(3))
graph save rd2.png,  replace
graph export rd2, replace

//graph  combine rd11.gph rd12.gph rd21.gph rd22.gph
//-2. Covaraites & Running variable-//
ds
foreach var of varlist tothsun8_nbr-blt40_yrs80occ_nbr{
	rdplot `var' hrs_82, c(28.5) p(3) plot
	graph export figure_`var'.png,replace
	rdrobust `var' hrs_82, c(28.5) p(3)
        }

//-3. Density of the running variable-//
rddensity hrs_82, c(28.5) p(1) plot //by default=median
```

**Sharp RD**

-  Treatment probability jumps at $c$ from 0 to 1

- $D_i$ is determined by $X_i$, i.e., $D_i=1(X_i\geq c)$

- We want to estimate $\tau_{SRD}=\mathbb{E}[Y_i(1)-Y_i(0)|X_i=c]$

- In practice, we implement 

  - For linear case: $Y_i=\alpha+\tau D_i+\beta(X_i-c)+\gamma (X_i-c)D_i+X_i+u_i$ s.t. $c-h<X_i<c+h$.

  - For polynomial case:

    $Y_i=\alpha+m(X_i)+\tau D_i+v_i$, where $m(X_i)$ include (1) polynomial of the running variable $X_i$, (2) fully interaction with the treatment effect $D_i$

  - If bandwith $h$ contains the range of $X_i$, then this is global polynomial; otherwise, local polynomial

**Fuzzy RD** 

- Treatment probability jumps at $c$
- We want to estimate $\tau_{FRD}=\frac{\lim_{x\rightarrow c-}\mathbb{E}[Y_i|X_i=x]-\lim_{x\rightarrow c+}\mathbb{E}[Y_i|X_i=x]}{\lim_{x\rightarrow c-}\mathbb{E}[D_i|X_i=x]-\lim_{x\rightarrow c+}\mathbb{E}[D_i|X_i=x]}$, 

- Further assumption:

  1. Monotonicity assumption: $D_i(x^{*})$ is non-increasing in $x^{*}$ at $x^{*}=c$  (An equivalent assumption in IV for '"no defiers"')
  2. Excludability: $X$ crossing the cutoff cannot impact receipt of treatment when these assumptions are made, which leads to  estimating $\mathbb{E}[Y_i(1)-Y_(0)|$ unit is  a complier and $X_i=c]$.

- In practice, we run the two stage regression and use the Wald estimator.

  - First stage: $Y_i=\pi_0+\pi_1Z_i+\pi_2(X_i-c)+\pi_3(X_i-c)Z_i+$Covariates$_i+u_i$, s.t., $c-h<X_i<c+h$.

  - Second stage: $D_i=\pi_0+\pi_1Z_i+\pi_2(X_i-c)+\pi_3(X_i-c)Z_i+$Covariates$_i+u_i$, s.t., $c-h<X_i<c+h$​.

  - $Z_i=1(x_I\geq C)$, $\hat{\tau}_{FRD}=\frac{\hat{\pi}_1}{\hat{\gamma}_1}$

  - If for polynomial, similar to the Sharp case. Same modification for the local and global.

  - ```stata
    //In rdrobust, there are three options for the kernel: triangular, epanechnikov, and uniform
    ////Local estimation
    //The variable in cluster must be a numerical
    //fuzzy(npl2000)  
    rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(1) fuzzy(npl2000) covs(tothsun8_nbr-hsdrop8_nbr  povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr) kernel(triangular) vce(cluster statefips) bwselect(mserd) 
    est store Q51
    
    ////Global estimation
    sum hrs_82
    local hvalueR=r(max)
    local hvalurL=abs(r(min))
    rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(3) fuzzy(npl2000) covs(tothsun8_nbr-hsdrop8_nbr povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr) h(`hvalueL'  `hvalueR') kernel(epanechnikov)  vce(cluster statefips)
    est store Q54
    ```

**Robustness Check**

1. Change the threshold for the running variable

   ```stata
   sum hrs_82
   local hrs_82max=r(max)
   local hrs_82min=abs(r(min))
   forvalues i=1(1)5{
   local jr=(`hrs_82max'-28.5)/(`i'+5)+28.5
   local jl=28.5-(28.5-`hrs_82min')/(`i'+5)
   rdrobust lnmdvalhs0_nbr hrs_82,c(`jr') p(3) fuzzy(npl2000) covs(tothsun8_nbr-hsdrop8_nbr povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr)  kernel(epanechnikov) vce(cluster statefips) bwselect(mserd) 
   estimates store jr`i'
   rdrobust lnmdvalhs0_nbr hrs_82,c(`jl') p(3) fuzzy(npl2000) covs(tothsun8_nbr-hsdrop8_nbr povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr)  kernel(epanechnikov) vce(cluster statefips) bwselect(mserd) 
   estimates store jl`i'
   }
   rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(3) fuzzy(npl2000) covs(tothsun8_nbr-hsdrop8_nbr  povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr) kernel(epanechnikov) vce(cluster statefips) bwselect(mserd) 
   estimates store jbaseline
   local vlist "jl1 jl2 jl3 jl4 jl5 jbaseline jr5 jr4 jr3 jr2 jr1"
   coefplot `vlist',  yline(0) drop(_cons) vertical 
   graph export placebotest1.png,replace
   ```

2. Donut test

   ```stata
   sum hrs_82
   local hrs_82max=r(max)
   forvalues i=1(1)5{
   local j=(`hrs_82max'-28.5)*0.05*`i'
   rdrobust lnmdvalhs0_nbr hrs_82 if abs(hrs_82-28.5)>`j', c(28.5) p(3) covs(tothsun8_nbr-hsdrop8_nbr povrat8_nbr-blt20_30yrs80occ_nbr) kernel(tri) vce(cluster statefips)
   estimates store obrob`i'
   }
   local vlist "obrob1 obrob2 obrob3 obrob4 obrob5"
   coefplot `vlist' , yline(0) drop(_cons) vertical
   graph export placebotest2.png,replace
   ```

3. Change the bandwidth

   ```stata
   rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(3) covs(tothsun8_nbr-hsdrop8_nbr povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr) kernel(tri)  vce(cluster statefips)    
   local h = e(h_l)   //the optimal bandwidth
   
   forvalues i=1(1)6{
   local hrobust=`h'*0.5*(`i')
   rdrobust lnmdvalhs0_nbr hrs_82, c(28.5) p(3) fuzzy(npl2000) h(`hrobust') covs(tothsun8_nbr-hsdrop8_nbr povrat8_nbr-blt20_30yrs80occ_nbr blt40_yrs80occ_nbr) kernel(tri)  vce(cluster statefips)   
   estimates store hrob`i'
   }
   local vlist "hrob1 hrob2 hrob3 hrob4 hrob5 hrob6"
   coefplot `vlist'  ,  yline(0) drop(_cons) vertical 
   graph export placebotest3.png,replace
   ```

   
