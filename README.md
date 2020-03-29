# SahulHumanSpread
Code and data files necessary for reproducing cellular-automaton model of human spread across Sahul

This code and these data reproduce the results in the following paper (currently in review):

BRADSHAW, CJA, K NORMAN, S ULM, AN WILLIAMS, C CLARKSON, J CHADŒUF, SC LIN, Z JACOBS, RG ROBERTS, MI BIRD, LS WEYRICH, S HABERLE, S O’CONNOR, B LLAMAS, TJ COHEN, T FRIEDRICH, P VETH, M LEAVESLEY, F SALTRÉ. In review. Rapid early peopling of Late Pleistocene Sahul. Proceedings of the National Academy of Sciences of the USA.

Contact:
Professor Corey J. A. Bradshaw
Matthew Flinders Fellow in Global Ecology
College of Science and Engineering
Flinders University
Adelaide, South Australia
e-mail: corey.bradshaw@flinders.edu.au
URL: http://GlobalEcologyFlinders.com

The R file 'SahulHumanSpreadGithub.SingleScenarioAvg' produces average scenario outputs over a set number of iterations. The user can choose the particulars of the scenario (e.g., underlying K~NPP relationship, entry time(s), entry point(s), spatial clustering, stochastic variances, minimum viable population thresholds, etc.)

The two zipped files should be decompressed and their files placed in the same directory as the R code.

The file 'matrixOperators.R' includes necessary functions and is sourced directly within the R code file.

GLOBAL SENSITIVITY ANALYSIS: Also included is the R file 'AusHumanSpreadGSAGithub.R' needed to do the global sensitivity analysis of the underlying parameters on the rate of saturation.

