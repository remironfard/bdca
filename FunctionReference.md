# Functions that work with the EEG struct #

## pop\_bdca\_mat2labels.m ##

Import the labels from a Matlab variable

```
EEG = pop_bdca_mat2labels(EEG,labels)

 input

    labels : can be the name of a variable, or a numeric vector
             with the labels in {0,1}

 output

    EEG.bdca.labels  :  the labels
 
```

## pop\_bdca\_event2labels.m ##

Import the labels from a EEG.epoch(...).event...

```
EEG = pop_bdca_event2labels(EEG,contrast,versus)

 input

    

 output

    EEG.bdca.labels  :  the labels
 
```

## pop\_bdca\_train.m ##

Train an EEG classifier. This function is a wrapper for the `bilinlogistregmultigp.m` program.

```
 synopsis:
 
    EEG = pop_bdca_train(EEG,R,sigma,lt,nut,ls,nus,supportframes,folds)
 
   input:
 
         EEG : The epoched EEGLAB dataset. The epoch labels must be present
               in EEG.bdca.labels. You can use pop_bdca_mat2labels(...) 
               to do that.
           R : The number of components.
       sigma : Regularization strength (std dev of parameters).
               low => more regularization, high  => less regularization.
               Typical values are between 0.001 and 1.
          lt : Temporal smoothing length scale in MILLISECONDS.
               Typical values are between 10 and 100.
         nut : Temporal smoothing shape parameter.
               Low => more ripple, High => less ripple.
               Typical value is nut = 2.5.
          ls : Spatial smoothness length scale. Only valid if the
               dataset has electrode locations. Set to [] to disable.
               Typical value = 0.1  (the unit is in EEGLAB relative head units).
         nus : Spatial smoothness ripple. Typical value = 100. Set to [] to disable.
   supportframes : a vector with the temporal indices (samples) 
                   within epoch to use for training. Set to [] to use the 
                   entire epoch for training. If set otherwise, then 
                   EEG.data(:,supportframes,:) will be used.
       folds : How many folds for crossvalidation. Set to 1 for training using the
               entire set. Set to e.g. 5 for 5-fold crossvalidation.
 
 
   output:
 
       if fold (see above) == 1
 
          EEG.bdca.cht.u = u : the spatial components.
          EEG.bdca.cht.t = t : the temporal components.
          EEG.bdca.cht.b = b : the logistic argument offset.
 
       if fold (see above) > 1
 
 
          EEG.bdca.cvstats.Az = A            : cross validated AUC.
          EEG.bdca.cvstats.loglik = loglik;  : cross validated log likelihood.

```

## pop\_bdca\_run.m ##

Classify EEG

```
Synopsis:

     EEG = pop_bdca_run(EEG)

Output:
 
     EEG.bdca.Ey    : Ey(n)  = E[y(n)]
     EEG.bdca.pot   : pot(n) = pi(Xn)
 
     if G is present (i.e. if you ran POP_BDCA_CONSTRICA)
 
        EEG.bdca.ica.S : Independent "source" across trials. 
```

## pop\_bdca\_constrica.m ##

Resolve the ambiguities by constrained independent component analysis.

```
synopsis:
 
    EEG = pop_bdca_constrica(EEG,alpha)
 
  input:
 
    alpha  :  Set it so that the algorithm converges smoothly.
              try e.g. alpha = 1000 
 
  output:
 
    EEG.bdca.ica.G      : ambiguity matrix
    EEG.bdca.ica.S      : independent component sources
    EEG.bdca.ica.alpha  : remember the alpha used
```

## pop\_bdca\_forward\_parafac.m ##

After running pop\_bdca\_constrica(...) use this this function to get the least-squares component scalp maps and time courses.

```
  synopsis:
 
    EEG = pop_bdca_forward_parafac(EEG)
 
  output:
 
    EEG.bdca.fwd.parafac  :  the 'Factors' output from PARAFAC ( >> help parafac )
 
  dependencies:
 
    NWAY toolbox
```



---


# Behind the scenes #

## bilinlogistregmultigp.m ##

Train a bilinear classifier with Gaussian Process regularization.

```
  synopsis:
 
     [w0,a,b,Ka,Kb] = bilinlogistregmultigp(x_train,y_train,R,sigw0)
     [w0,a,b,Ka,Kb] = bilinlogistregmultigp(x_train,y_train,R,sigw0,gpa,gpb,spaceunits)
 
  inputs:
 
      x_train  : D-by-T-by-N data
      y_train  : 1-by-N, binary labels in {0,1}
            R  : number of components
        sigw0  : prior stddev of w0
      gpa,gpb  : prior covariance parameters for a and b resp.
                 : if left out, defaults to sigw0.
                 : otherwise, a scalar gives the prior stddev.
                 : otherwise, a vector [sigma l] defines a GP with SE cov. function.
                 : otherwise, a vector [sigma l nu] defines a GP with Matern cov. function
   spaceunits  : Spatial coords in euclidian basis [x y z; x y z; etc].
                 Set it to [] if you don't want to declare spatial smoothness, or if
                 you don't know electrode coordinates.
 
    outputs:
 
         w0,a,b : you know.
          Ka,Kb : covariancematrices for a and b respectively.  
 
  Dependencies: OGP, IMMOPTIBOX
```

## bilinlogistregmultigp\_run.m ##

Classify

```
[Ey,pot] = bilinlogistregmultigp_run(x,w0,a,b)
```

## auc.m ##

Compute the area under the ROC curve.

```
A = auc(L,fi)
 
    L  : a vector with the true labels in {0,1}
    fi : a vector with the "soft" predictions
```

## bilinlogistregmultigp\_bayes.m ##

Bayesian model selection / Marginal Likelihood. I haven't tested this function yet! - The cholesky might fail, the Hessian is not necessarily positive definite. But, BIC should work, no matter what.

In fact, a procedure for model tuning could be this:


  1. Select (and fix) the GP hyperparameters using Laplace with R=1. BIC is no good for this purpose, since it does not depend on the prior.
  1. Then increase R. Laplace is more likely to fail with R>1 hence use BIC in this step.

```
  [BIC,Laplace] = bilinlogistregmultigp_bayes(w0,a,b,X3,y,sigw0,Kl,Kr)
 
    or to avoid the Cholesky of the Hessian:
 
  BIC = bilinlogistregmultigp_bayes(w0,a,b,X3,y,sigw0,Kl,Kr)
```