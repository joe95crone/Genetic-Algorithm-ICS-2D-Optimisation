__Genetic Algorithm ICS 2D Optimisation__
2D Non-round beam collimated flux vs rms bandwidth optimisation using the SPEA2 algorithm on the PISA platform. 

This code is a minimisation algorithm designed by Balsa Terzic and Joe Crone which uses the SPEA2 minimisation algorithm to maximise the collimated flux (-Fcol) within an *RMS* bandwidth (BW). Variable sets are the IP beta-functions and the collimation angle. This fully takes into account recoil, angular crossing and hourglass effects within the linear regime. It uses my personal derived analytical collimated flux equation, which only neglects energy spread of the bunch and pulse, and uses the *RMS* bandwidth derivation by Ranjan et al. 

Currently this can only be used to create tuning curves, customization is ongoing for a penalty function version in which a single user specified BW point can be optimised.

Customization to include more optimisation variables, for example for a longitudinal optimisation and transverse optimisation simultaneously or square/rectangular collimation with collimation angles in each plane.
