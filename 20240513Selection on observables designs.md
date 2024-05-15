## Selection on observables

### Key assumptions 

1. CI (the conditional-independence assumption) $Y_i(0), Y_i(1) \perp D_i|X_i$ : the potential outcomes are independent of the treatment conditional on observable covariates
2. SUTVA (stable unit treatment value assumption/ no interference/ the i.i.d. sampling assumption): unit $i$'s potential outcomes are unaffected by whether unit $j (j\not=i)$​ is treated or untreated.

### Matching

- Ideal of matching: to compare treated units $(D_i=1)$ to control units $(D_i=0)$ that have similar values of $X_i$, so that we can estimate the treatment effect

$$
\tau(x)=\mathbb{E}[Y_i(1)-Y_i(0)|X_i=x],
$$

because given the selection on observables assumption, we know that $D_i$ is as good as randomly assigned after conditioning on $X_i$.

- Assume there are $N_T$ treated units and $N_C$ control units. Define $N_T$ sets of weights, with $N_C$ weights in each set:
  $$
  w_i(j)=(i=1,...,N_T, j=1,...,N_C)
  $$
  For each set of weights, let $\sum_jw_i(j)=1$. The generic matching estimator is
  $$
  \hat{\tau}_M=\frac{1}{N_T}\sum_{i\in \{D=1\}}[y_i-\sum_{j\in \{D=0\}}w_i(j)y_j].
  $$

  - If $w_i(j)=\frac{1}{N_C}$, then $\hat{\tau}_M=\frac{1}{N_C}$
  - Nearest-neighbor matching: set $w_i(j)$ as the distance between $X_i$ and $X_j$, with more weight on the closer pairs and less weight on the farther pairs.

- But in the matching, we face the curse of dimensionality: the more variables we have in $X$, the less likely we can find a comparison control unit lying close to any given treatment unit.

### Propensity Score

Useful manual: [teffects](https://www.stata.com/manuals/teteffectsintro.pdf)

- Further assumption:  [the overlap assumption (also called the "strongly ignorable" or "common support") holds](https://www.urban.org/research/data-methods/data-analysis/quantitative-data-analysis/impact-analysis/quasi-experimental-methods/propensity-score-matching): $0<P(D_i=1|X_i)<1$ such that for each group divided by $X$​, we have treated units and untreated units. 

  - [Test for overlap in Stata](https://www.stata.com/manuals14/teteffectsoverlap.pdf)

    ```stata
    use http://www.stata-press.com/data/r14/cattaneo2
    teffects ipw (bweight) (mbsmoke mmarried c.mage##c.mage fbaby medu, probit)
    teffects overlap
    ```

  - Under the "strongly ignorable" treatment assumption, it is sufficient to condition simply on $p(X_i)=\mathbb{E}[D_i|X_i]$, known as the **propensity score**: if $(Y_i(0), Y_i(1))\perp D_i |X_i$, then $(Y_i(0), Y_i(1))\perp D_i | p(X_i)$.

- How to get the propensity score? 
  - [Code example](https://www.lianxh.cn/details/173.html)

```stata
logit treat x //use the logit model 
//or you can use probit treat x
est store model1 //save the model
predict p //get the propensity score
label var p propensity_score
```

- What can we do with the propensity score? 

  1. Regression adjustment on the propensity score (include the propensity score as a regressor)

     - If we believe that the treatment effects are homogeneous

       $Y_i=\alpha+\delta D_i+\beta\hat{p}_i+u_i$,​

       the estimated treatment effect for any given $X_i$ is $\delta$.

     ```stata 
     reg y treat p, vce(cluster stresfip)
     est store m33
     ```

     - If we believe that the treatment effects are heterogeneous

       $Y_i=\alpha+\delta_1D_i+\delta_2D_i\hat{p}_i+\beta\hat{p}_i+u_i$

       the estimated treatment effect for any given $X_i$ is $\delta_1+\delta_2 \hat{p}_i$, and the average treatment effect is $\delta_1+\delta_2 \bar{p}$

       ```stata 
       reg y treat p treat#c.p
       //or reg y i.treat##c.p
       ```

       Note: how to incorporate [the interaction term](https://bbs.pinggu.org/thread-8628238-1-1.html)

  2. Blocking

  - Idea: divide the range of the propensity score into $K$ blocks. within each block $k$, compute $\hat{\tau}_k$. Finally, combine all $K$ treatment effect estimates as follows: $$\hat{\tau}=\sum_{k=1}^K \hat{\tau}_k \frac{N_k}{N}=\sum_{k=1}^K\hat{\tau}_k\frac{N_{1k}+N_{0k}}{N}$$

  - We can see whether the overlap assumption holds through blocking.

  - When using blocking, one should discard all control units with $\hat{p}(X_i)$ less than the minimum $\hat{p}(X_i)$ in the treated group and all treated units with a $\hat{p}(X_i)$ greater than the maximum $\hat{p}(X_i)$ in the control group.  e.g.

    ```stata
    by smoke: egen min_p0=min(p0)
    egen cutoff_min_p0=max(min_p0)
    by smoke: egen max_p0=max(p0)
    egen cutoff_max_p0=min(max_p0)
    drop if p0<cutoff_min_p0 | p0>cutoff_max_p0
    ```

  3. Reweighting

     What we apply in practice is the inverse-probability weighting ([IPW](https://www.stata.com/manuals/teteffectsipw.pdf#teteffectsipw)).
     $$
     \hat{\tau}=(\frac{\sum_{i=1}^N\frac{D_iY_i}{\hat{p}(X_i)}}{\sum_{i=1}^N\frac{D_i}{\hat{p}(X_i)}})-(\frac{\sum_{i=1}^N\frac{(1-D_i)Y_i}{1-\hat{p}(X_i)}}{\sum_{i=1}^N\frac{1-D_i}{1-\hat{p}(X_i)}})
     $$

     ```stata
     //e.g. We estimate the effect of smoking (treatment) on birth weight (y) by using a probit model to predict the mother's smoking behavior as a function of X
     teffects ipw (bweight) (mbsmoke mmarried mage prenatal1 fbaby, probit)
     ```

  4. Matching

     The famous PSM (propensity-score matching): matches on the estimated predicted probabilities of treatment. PSM does not require bias correction.

     ```stata
     //Here we model the propensity score using a probit model, incorporating marital status...
     teffects psmatch (bweight) (mbsmoke mmarried mage prenatal1 fbaby, probit)
     ```
