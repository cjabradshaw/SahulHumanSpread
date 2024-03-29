# Sahul Human Spread
<img align="right" src="ausdots.png" alt="Australia cellular" width="200" style="margin-top: 20px">

Code and data files necessary for reproducing cellular-automaton model of human spread across Sahul <a href="http://doi.org/10.5281/zenodo.4453767">
  
  <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4453767.svg"></a>

This code and these data reproduce the results in the following paper:

<a href="http://www.flinders.edu.au/people/corey.bradshaw">BRADSHAW, CJA<a/>, <a href="https://www.linkedin.com/in/kasih-norman-7a122a110/?originalSubdomain=au">K NORMAN</a>, <a href="https://research.jcu.edu.au/portfolio/sean.ulm">S ULM</a>, <a href="https://www.emmconsulting.com.au/about/leadership-team/dr-alan-william-2/">AN WILLIAMS</a>, <a href="http://researchers.uq.edu.au/researcher/742">C CLARKSON</a>, <a href="https://www.researchgate.net/profile/Joel_Chadoeuf">J CHADŒUF</a>, <a href="https://scholars.uow.edu.au/display/sam_lin">SC LIN</a>, <a href="https://scholars.uow.edu.au/display/zenobia_jacobs">Z JACOBS</a>, <a href="https://scholars.uow.edu.au/display/richard_roberts">RG ROBERTS</a>, <a href="https://research.jcu.edu.au/portfolio/michael.bird">MI BIRD</a>, <a href="https://anth.la.psu.edu/people/lsw132">LS WEYRICH</a>, <a href="https://researchers.anu.edu.au/researchers/haberle-sg">S HABERLE</a>, <a href="https://researchers.anu.edu.au/researchers/o-connor-sl">S O’CONNOR</a>, <a href="https://www.adelaide.edu.au/directory/bastien.llamas">B LLAMAS</a>, <a href="https://scholars.uow.edu.au/display/tim_cohen">TJ COHEN</a>, <a href="http://iprc.soest.hawaii.edu/users/tobiasf/">T FRIEDRICH</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-veth">P VETH</a>, <a href="https://www.upng.ac.pg/index.php/shss-staff-division/aas-shss-contact/144-dr-matthew-leavesley">M LEAVESLEY</a>, <a href="http://www.flinders.edu.au/people/frederik.saltre">F SALTRÉ</a>. 2021. <a href="http://doi.org/10.1038/s41467-021-21551-3">Stochastic models support rapid peopling of Late Pleistocene Sahul</a> <em>Nature Communications</em> 12: 2440. doi:10.1038/s41467-021-21551-3

## Abstract
The peopling of Sahul (the combined continent of Australia and New Guinea) represents the earliest continental migration and settlement event of solely anatomically modern humans, but its patterns and ecological drivers remain largely conceptual in the current literature. We present an advanced stochastic-ecological model to test the relative support for scenarios describing where and when the first humans entered Sahul, and their most probable routes of early settlement. The model supports a dominant entry via the northwest Sahul Shelf first, potentially followed by a second entry through New Guinea, with initial entry most consistent with 50,000 or 75,000 years ago based on comparison with bias-corrected archaeological map layers. The model’s emergent properties predict that peopling of the entire continent occurred rapidly across all ecological environments within 156–208 human generations (4368–5599 years) and at a plausible rate of 0.71–0.92 km year<sup>−1</sup>. More broadly, our methods and approaches can readily inform other global migration debates, with results supporting an exit of anatomically modern humans from Africa 63,000–90,000 years ago, and the peopling of Eurasia in as little as 12,000–15,000 years via inland routes.
  
Contact:
Professor Corey J. A. Bradshaw, Matthew Flinders Professor of Global Ecology

College of Science and Engineering, Flinders University

Adelaide, South Australia

  e-mail: <a href="mailto:corey.bradshaw@flinders.edu.au">corey.bradshaw@flinders.edu.au</a>; URL: http://GlobalEcologyFlinders.com

The R file <a href="https://github.com/cjabradshaw/SahulHumanSpread/blob/master/SahulHumanSpreadGithub.SingleScenarioAvg.R"><code>SahulHumanSpreadGithub.SingleScenarioAvg.R</code></a> produces average scenario outputs over a set number of iterations. The user can choose the particulars of the scenario (e.g., underlying <em>K</em>~NPP relationship, entry time(s), entry point(s), spatial clustering, stochastic variances, minimum viable population thresholds, etc.)

The two zipped files should be decompressed and their files placed in the same directory as the R code.

The file <a href="https://github.com/cjabradshaw/SahulHumanSpread/blob/master/matrixOperators.r"><code>matrixOperators.R</code></a> includes necessary functions and is sourced directly within the R code file.

The file '<a href="https://github.com/cjabradshaw/SahulHumanSpread/blob/master/Archaeology%20sites%20%26%20dates%20used%20for%20comparison%20layers.xlsx">Archaeology sites & dates used for comparison layers.xlsx</a>' is an Excel file listing all the archaeological specimen dates, type of material, dating metthod, dating technique, and quality rating used in constructing the archaeological comparison layers. 

NOTE: The code to reproduce these spatial layers can be sourced from <a href="https://github.com/FredSaltre/SEOZ_megafauna_extirpation">this Github repository</a> (associated with the paper: Saltré, F, J Chadoeuf, KJ Peters, MC McDowell, T Friedrich, A Timmermann, S Ulm, CJA Bradshaw. 2019. <a href="http://doi.org/10.1038/s41467-019-13277-0">Climate-human interaction associated with southeast Australian megafauna-extinction patterns</a>. <em>Nature Communications</em> 10: 5311. doi:10.1038/s41467-019-13277-0)

<strong>GLOBAL SENSITIVITY ANALYSIS</strong>: Also included is the R file <a href="https://github.com/cjabradshaw/SahulHumanSpread/blob/master/AusHumanSpreadGSAGithub.R"><code>AusHumanSpreadGSAGithub.R</code></a> needed to do the global sensitivity analysis of the underlying parameters on the rate of saturation.

