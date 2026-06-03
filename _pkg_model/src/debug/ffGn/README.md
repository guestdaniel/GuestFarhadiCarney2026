# Summary
At least two (possibly more?) versions of the `ffGn.m` MATLAB function to synthesize fractional Gaussian noise have existed at different points in time to be used in different versions of the auditory-nerve model.

# Sources
- `ffGn_2014.m`: [Zilany et al. (2014) model code](https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/codes/Zilany-2014-Code-and-paper.zip)
- `ffGn_urear_2020b.m`: [UR Ear 2020b code](https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/codes/UR_EAR_2020b.zip)

# Differences
Both versions of the code have the same arguments, (`N`, `tdres`, `Hinput`, `noiseType`, `mu`, and `sigma`), and for the same values of those arguments, behave in the exact same way.
However, at least in the 2014 model, `ffGn` is called via a C-to-MATLAB callback with only the first five positional arguments (`N`, `tdres`, `Hinput`, `noiseType`, `mu`), leaving the last to be set based on argument-handling logic inside of `ffGn.m`.
In the 2014 code, values of `0.1`, `4.0`, and `100.0` for passed in for `mu` (based on the `fiberType` switch) to produce LSR, MSR, and HSR fibers, respectively.
The 2014 version of the code sets the value of `sigma` according to the following logic based on the value of `mu`:
```
if mu<0.5
    sigma = 3;%5  
else
    if mu<18
        sigma = 30;%50   % 7 when added after powerlaw
    else
        sigma = 200;  % 40 when added after powerlaw        
    end
end
```

In contrast, the 2018 version of the code sets the value of `sigma` differently, to compensate for other changes made in the 2014-to-2018 model update:
```
if mu<0.2
    sigma = 1;%5  
else
    if mu<20
        sigma = 10;
    else
        sigma = mu/2;
    end
end
```
In fact, the only time `mu` is even used is to set `sigma`!
The final result is that values of `3.0`, `30.0`, and `200.0` are used as values for `sigma` for LSR, MSR, and HSR fibers, respectively, in the 2014 version of the code, whereas values of `1.0`, `10.0`, and `50.0` are used as values in the 2018 version of the code, assuming that the correct `mu` values are passed into the function.
In summary, if the 2018 ffGn code is used with the 2014 model, the result is an underestimate of the randomness of mean rates from simulation from simulation, both for spontaneous and (to a lesser extent) driven rates.

The efferent model also makes calls to `ffGn.m`, and it turns out I've been using the wrong `ffGn.m` with the efferent model (and bundled the wrong one in the last release I sent to everyone --- sorry!).
However, I pass a `sigma` value when I call `ffGn` (i.e., I pass all 6 arguments) instead of relying on internal logic to set `sigma`, so users of the efferent model should be unaffected by which version of the code is on their path.
That said, investigating this helped me discover I was passing in erroneous values for `mu` and `sigma` to the efferent model, which also resulted in ffGn with insufficient variance to match the 2014 model --- a different, but related, bug!
This bug has been fixed in the attached version, as has the `ffGn.m` function with the updated random number generation from last week, and the efferent model should now be getting the "correct" ffGn regardless of what version of `ffGn.m` is called.

