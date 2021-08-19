__Genetic Algorithm ICS 2D Optimisation__

2D Non-round beam collimated flux vs rms bandwidth optimisation using the SPEA2 algorithm on the PISA platform. 

This code is a minimisation algorithm designed by Balsa Terzic and Joe Crone which uses the SPEA2 minimisation algorithm to maximise the collimated flux (-Fcol) within an *RMS* bandwidth (BW). Variable sets are the IP beta-functions and the collimation angle. This fully takes into account recoil, angular crossing and hourglass effects within the linear regime. It uses my personal derived analytical collimated flux equation, which only neglects energy spread of the bunch and pulse, and uses the *RMS* bandwidth derivation by Ranjan et al. 

Customization to include more optimisation variables, for example for a longitudinal optimisation and transverse optimisation simultaneously or square/rectangular collimation with collimation angles in each plane.

A tuning curve version (TuningCurve_NO_DAT) has been produced, which produces tuning curves within a certain bandwidth range specified (rarther roughly) by the collmation angle variable.

A single point version (SinglePoint_NO_DAT) which uses a penalty function to focus on a target bandwidth has been produced. The maximised flux point can either be found by taking the minimum value of the simulation or plotting the y intercepts of the data against the penalty function.

The tuning curve version has also been modified to produce results for the case where the laser pulse and electron bunch are crabbed (CC_TuningCurve_NO_DAT).  


