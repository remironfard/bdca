# Before reading this tutorial #

The tutorial uses an example dataset. You can download it [here](http://code.google.com/p/bdca/downloads/list). When you have the example dataset proceed as follows...


---


# Tutorial for using BDCA with EEGLAB #

## Step 1. Epoch your data ##

This is standard EEGLAB stuff. If you load the dataset in EEGLAB and scroll through it, you will see a bunch of events. We'll use those to epoch the dataset. You can do it from the EEGLAB menus, or simply run the following from the command line:
```
EEG = pop_epoch( EEG, {  '80'  '160'  }, [-0.1         0.5], 'newname', 'BDCA tutorial dataset epochs', 'epochinfo', 'yes');
```
So, now we have epochs available.

_EEGLAB BUG: Note that if you set 0 instead of -0.1 above then EEGLAB will make a sometimes/somewhat useless event structure. Therefore, always include some samples preceeding the event, e.g. -0.1 instead of 0._

## Step 2. Make BDCA know which labels to use for what ##

Ok, so we want to train a classifier which discriminates between epochs of label 160 and epochs of label 80. Now we use the BDCA toolbox:
```
EEG = pop_bdca_event2labels(EEG,160,80);
```
This tells the classifier that you want to discriminate epochs of type 160 against epochs of type 80. What happened is that a vector is now to be found in `EEG.bdca.labels` !  take a look at it if you want:
```
figure, plot(EEG.bdca.labels);
```

Note: if you wanted to discriminate, say 160 AND 170 versus 80 then you could do e.g.
```
EEG = pop_bdca_event2labels(EEG,{160 170},80);
```
... but we don't need that for this tutorial.

You can manually set the labels by changing `EEG.bdca.labels`. The function `pop_bdca_mat2labels` is meant to help you with that.

## Step 3. Tuning sigma / scaling the data ##
We're gonna tune `sigma` using 5-fold cross validation. Our initial guess is `sigma=0.1`, and also let's set `R=1` because it is much faster than `R>1` and it is sufficient for tuning `sigma`:
```
R=1;
sigma=0.1;
EEG = pop_bdca_train(EEG,R,sigma,[],[],[],[],[],5);
```
when it is done, you can see the performance ROC area in the variable `EEG.bdca.cvstats.Az`. The cross validated log-likelihood is also there in `EEG.bdca.cvstats.loglik`. So you can tune your model looking at either one.

So, now you could imagine tuning this by hand, or making a for loop and choosing `sigma` based on the Az or loglik. Note that sigma should be sweeped on a logarithmic scale, say
```
happy = [];
for sigma=logspace(-2,0,30)
  EEG = pop_bdca_train(EEG,R,sigma,[],[],[],[],[],5);
  happy = [happy; sigma,EEG.bdca.cvstats.Az];
end
```
and then choosing the `sigma` which has the highest `Az` in the result `happy`. Or, like I said, use the log-likelihood.

Due to numerical issues, you might wanna scale your data so that the optimal sigma is roughly 1. Try for instance to scale your data by the reciprocal of the average channel standard deviation.

## Step 4. Tuning temporal smoothness ##

Now the time has come to start tuning the temporal smoothness. I recommend setting `nut=2.5` but you might find a different value to love. Anyway, tune `lt`...
```
lt=30;  % this means that we expect 30 milliseconds smoothness
nut=2.5;
EEG = pop_bdca_train(EEG,R,sigma,lt,nut,[],[],[],5);
```
and, again, look at `EEG.bdca.cvstats.Az` or `loglik`. As before with `sigma`, make a for-loop for finding the best `lt`.

## Step 5. Tuning spatial smoothness ##

Tuning the spatial smoothness **requires that you have channel locations for all channels and that no two channels have the same coordinates**. (I'm just saying so that you don't get tempted and simply copy locations from other electrodes to channels that have no coordinates).

```
ls=0.1;  % this is relative to the EEGLAB head radius which I believe is 0.5
nus=100; % 100 means: no ripple. Decrease this number to allow more ripply topographys
EEG = pop_bdca_train(EEG,R,sigma,lt,nut,ls,nus,[],5);
```
So, tune `ls` and `nus` using `cvstats` as before.

## Step 6. Increasing the number of components ##

Finally, try increasing `R` to see the effect of having more components. Note that it will start to go pretty slow.

## Step 7.  ... much more to come ... ##



---

# Tutorial for using BDCA without EEGLAB #

## Tips on tuning the parameters with bilinlogistregmulti ##

First, set R to 1
```
>> R = 1
```
Set sigw0 to some high number, like 5
```
>> sigw0 = 5;
```
Then, tune sigma for space and time
```
>> gpa = [sigma];
>> gpb = [sigma];
>> spaceunits = [];
>> [w0,a,b] = bilinlogistregmultigp(X,y,R,sigw0,gpa,gpb,spaceunits);
```
repeat the above step until you find a good sigma (use cross validation!). Then, start tuning the temporal smoothness. I recommend setting nut to 2.5, and then try different values for lt.
```
>> nut = 2.5;
>> gpb = [sigma lt nut];
>> [w0,a,b] = bilinlogistregmultigp(X,y,R,sigw0,gpa,gpb,spaceunits);
```
repeat the above step until you find a good lt. If you have channel locations, then you can tune the spatial smoothness. First, you have to set spaceunits non-empty. Then start tuning nus and ls.
```
>> gpa = [sigma ls nus];
```
Finally, you can start to increase R in increments of 1.