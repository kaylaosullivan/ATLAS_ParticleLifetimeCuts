# ATLAS_ParticleLifetimeCuts
Automatically apply cuts on available parameters to determine which cuts reduce the background-to-signal ratio the most.

Code is currently set up to perform cuts and load Lambda and Kaon particle data. This can easily be modified.


Use performCuts(...) to try multiple values for cuts. For each cut, the script will:
- Refit the data with two exponentials
- Calculate the lifetime of the particle based on the new fit
- Calculate the background-to-signal ratio that we're trying to reduce according to the fit
- Output all information relating to the cut, background:signal ratios, lifetimes, fit parameters, etc to a csv file

Use errorCuts(...) to do as mentioned above, with the addition of testing the upper and lower bounds of each cut to determine its effect on lifetime.
