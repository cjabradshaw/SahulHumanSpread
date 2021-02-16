# SahulHumanSpread
Code and data files necessary for reproducing cellular-automaton model of human spread across Sahul <a href="http://doi.org/10.5281/zenodo.4453767">
  
  <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4453767.svg"></a>

This code and these data reproduce the results in the following paper:

<a href="http://www.flinders.edu.au/people/corey.bradshaw">BRADSHAW, CJA<a/>, <a href="https://www.linkedin.com/in/kasih-norman-7a122a110/?originalSubdomain=au">K NORMAN</a>, <a href="https://research.jcu.edu.au/portfolio/sean.ulm">S ULM</a>, <a href="https://www.emmconsulting.com.au/about/leadership-team/dr-alan-william-2/">AN WILLIAMS</a>, <a href="http://researchers.uq.edu.au/researcher/742">C CLARKSON</a>, <a href="https://www.researchgate.net/profile/Joel_Chadoeuf">J CHADŒUF</a>, <a href="https://scholars.uow.edu.au/display/sam_lin">SC LIN</a>, <a href="https://scholars.uow.edu.au/display/zenobia_jacobs">Z JACOBS</a>, <a href="https://scholars.uow.edu.au/display/richard_roberts">RG ROBERTS</a>, <a href="https://research.jcu.edu.au/portfolio/michael.bird">MI BIRD</a>, <a href="https://anth.la.psu.edu/people/lsw132">LS WEYRICH</a>, <a href="https://researchers.anu.edu.au/researchers/haberle-sg">S HABERLE</a>, <a href="https://researchers.anu.edu.au/researchers/o-connor-sl">S O’CONNOR</a>, <a href="https://www.adelaide.edu.au/directory/bastien.llamas">B LLAMAS</a>, <a href="https://scholars.uow.edu.au/display/tim_cohen">TJ COHEN</a>, <a href="http://iprc.soest.hawaii.edu/users/tobiasf/">T FRIEDRICH</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-veth">P VETH</a>, <a href="https://www.upng.ac.pg/index.php/shss-staff-division/aas-shss-contact/144-dr-matthew-leavesley">M LEAVESLEY</a>, <a href="http://www.flinders.edu.au/people/frederik.saltre">F SALTRÉ</a>. 2021. <a href="http://doi.org/10.1038/s41467-021-21551-3">Stochastic models support rapid peopling of Late Pleistocene Sahul</a> <em>Nature Communications</em>. doi:10.1038/s41467-021-21551-3

Contact:
Professor Corey J. A. Bradshaw, Matthew Flinders Fellow in Global Ecology

College of Science and Engineering, Flinders University

Adelaide, South Australia

e-mail: corey.bradshaw@flinders.edu.au; URL: http://GlobalEcologyFlinders.com

The R file 'SahulHumanSpreadGithub.SingleScenarioAvg' produces average scenario outputs over a set number of iterations. The user can choose the particulars of the scenario (e.g., underlying K~NPP relationship, entry time(s), entry point(s), spatial clustering, stochastic variances, minimum viable population thresholds, etc.)

The two zipped files should be decompressed and their files placed in the same directory as the R code.

The file 'matrixOperators.R' includes necessary functions and is sourced directly within the R code file.

The file 'Archaeology sites & dates used for comparison layers.xlsx' is an Excel file listing all the archaeological specimen dates, type of material, dating metthod, dating technique, and quality rating used in constructing the archaeological comparison layers. NOTE: The code to reproduce these spatial layers can be sourced from https://github.com/FredSaltre/SEOZ_megafauna_extirpation (associated with the paper: Saltré, F, J Chadoeuf, KJ Peters, MC McDowell, T Friedrich, A Timmermann, S Ulm, CJA Bradshaw. 2019. Climate-human interaction associated with southeast Australian megafauna-extinction patterns. Nature Communications 10: 5311. doi:10.1038/s41467-019-13277-0)

GLOBAL SENSITIVITY ANALYSIS: Also included is the R file 'AusHumanSpreadGSAGithub.R' needed to do the global sensitivity analysis of the underlying parameters on the rate of saturation.

