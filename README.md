---
author: Nicola F. Müller
level: Intermediate
title: MASCOT Skyline v3.0.5 Tutorial
subtitle: Parameter and State inference using the approximate structured coalescent
beastversion: 2.7.x
tracerversion: 1.7.2
---


# Background

The structured coalescent models how lineages coalesce within and migrate between sub-populations using coalescent and migration rates, which can be related to epidemiologically more meaningful parameters, such as the prevalence and transmission rates {% cite volz2012complex --file MascotSkyline-Tutorial/master-refs %}. Structured coalescent methods largely assume that the rates of coalescence and migration are constant over time, though deterministic approaches to model parametric dynamics from compartmental models exist {% cite volz2018bayesian --file MascotSkyline-Tutorial/master-refs %}. Here, we introduce a phylodynamic framework to infer non-parametric effective population size (Ne) dynamics under the marginal approximation of the structured coalescent (MASCOT). The effective population sizes are estimated at predefined points in time, between which we assume exponential growth dynamics\citep{volz2018modeling}. As such, we allow the Ne's to continuously change over time instead of assuming piecewise constant dynamics, as is typically used in skyline approaches (for example~\citep{gill2013improving}). We use a Gaussian Markov Random Field (GMRF), as in~\citep{gill2013improving} for unstructured populations, to model the temporal correlation Ne's. We then estimate Ne trajectories for each sub-population in the model using Markov chain Monte Carlo (MCMC) by using MCMC operations that learn the correlations structure between the different parameters~\cite{baele2017adaptive}. We first show, using simulations, that we can retrieve non-parametric population dynamics and migration rates of different sub-populations from phylogenetic trees. We then show how accounting for population structure improves the inference of population dynamics and vice versa. Lastly, we compare the ancestral state reconstruction and inference results of migration rates between MASCOT-Skyline and DTA~\citep{lemey2009bayesian} using a dataset of SARS-CoV-2 sequences and Susceptible-Infected-Recovered (SIR) simulations. We implemented MASCOT-Skyline as part of the BEAST2 package MASCOT~\citep{muller2018mascot}, an addon to the Bayesian phylogenetics software platform BEAST2~\citep{bouckaert2019beast}

----

# Programs used in this Exercise

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 ([http://www.beast2.org](http://www.beast2.org)) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees. This tutorial is written for BEAST v{{ page.beastversion }} {% cite BEAST2book2014 --file Mascot-Tutorial/master-refs %}.


### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### TreeAnnotator

TreeAnnotator is used to summarise the posterior sample of trees to produce a maximum clade credibility tree. It can also be used to summarise and visualise the posterior estimates of other tree parameters (e.g. node height).

TreeAnnotator is provided as a part of the BEAST2 package so you do not need to install it separately.


### Tracer

Tracer ([http://tree.bio.ed.ac.uk/software/tracer](http://tree.bio.ed.ac.uk/software/tracer)) is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v{{ page.tracerversion }}.


### IcyTree.org

IcyTree.org is a web-based tree viewer that, in addition to trees, also allows the user to visualize networks in the extended newick format {% cite vaughan2017icytree --file Mascot-Tutorial/master-refs %}.

----

# Practical: Parameter and State inference using the approximate structured coalescent

In this tutorial, we will estimate migration rates, nonparametric (skygrid) population sizes and the locations of internal nodes using the marginal approximation of the structured coalescent implemented in BEAST2, MASCOT {% cite mueller2017mascot --file Mascot-Tutorial/master-refs %}.

The aim is to:

-  Learn how to set up a MASCOT-Skyline analysis
-  Learn how to infer structure from trees with sampling location
-  Learn how to estimate time-varying effective population sizes and migration rates
-  Get to know how to choose the set-up of such an analysis
-  Learn how to read the output of a MASCOT-Skyline analysis

## Setting up an analysis in BEAUti

### Download MASCOT

First, we have to download the package MASCOT (we need at least version 3.0.5) using the BEAUTi package manager. Go to _File >> Manage Packages_ and download the package MASCOT.

<figure>
	<a id="fig:example1"></a>
	<img style="width:50%;" src="figures/MascotDownload.png" alt="">
	<figcaption>Figure 1: Download the MASCOT package.</figcaption>
</figure>

MASCOT will only be available in BEAUti once you close and restart the program.

### Loading the Influenza A/H3N2 Sequences (Partitions)

The sequence alignment is in the file [H3N2.nexus](http://github.com/nicfel/MascotSkyline-Tutorial/raw/master/data/sequences.nexus).
Right-click on this link and save it to a folder on your computer.
Once downloaded, this file can either be drag-and-dropped into BEAUti or added by using BEAUti's menu system via _File >> Import Alignment_.
Once the sequences are added, we need to specify the sampling dates and locations.

### Get the sampling times (Tip Dates)

Open the "Tip Dates" panel and then select the "Use tip dates" checkbox.

The sampling times are encoded in the sequence names.  We can tell BEAUti to use these by clicking the _Auto-configure_ button. The sampling times are in the isolate name after the last vertical bar "|" in the sequence name. To extract these times, select "use everything", select "after last" and enter "|" (without the quotes) in the text box immediately to the right. The setup should look as shown in the figure below.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/TipDates.png" alt="">
	<figcaption>Figure 2: Guess sampling times.</figcaption>
</figure>

Clicking "OK" should now populate the table with the sample times extracted from the sequence names: the column **Date** should now have values between 2016 and 2015 and the column **Height** should have values from 0 to 2. The heights denote the time difference from a sequence to the most recently sampled sequence. If everything is specified correctly, the sequence with Height 0.0 should have a date of 2016.77.

### Specify the Site Model (Site Model)

Next, we have to specify the site model. To do this, choose the "Site Model" tab. For Influenza Hemagluttanin sequences as we have here, HKY is the most commonly used model of nucleotide evolution. This model allows for differences in transversion and transition rates, meaning that changes between bases that are chemically more closely related (transitions) are allowed to have a different rate to changes between bases that are chemically more distinct (transversions).
Additionally, we should allow for different rate categories for different sires in the alignment.
This can be done by setting the _Gamma Category Count_ to 4, which is just a value that has typically been used. Make sure that the estimate is checked next to the shape parameter. To reduce the number of parameters we have to estimate, we can set Frequencies to Empirical.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/SiteModel.png" alt="">
	<figcaption>Figure 4: Set the site model.</figcaption>
</figure>


### Set the clock model (Clock Model)

For rapidly evolving viruses, the assumption of a strict molecular clock is often made, meaning that the molecular clock is the same on each branch of the phylogeny. We will leave everything to the default here

### Specify the priors (Priors)

We first have to choose the tree prior, which in this case is MASCOT.
To do so, search the drop-down menu next to `Tree.t:sequences` and choose MASCOT.
To see more options, we have to expand the MASCOT tree prior by clicking the arrow to the left of the label.
By default, the rate dynamics for this setting are `Constant`, which means that effective population sizes and migration rates are assumed to be constant through time.
To use skyline dynamics, we have to choose `Skyline` from the drop-down menu next to `Dynamics`.
<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/Skyline.png" alt="">
	<figcaption>Figure 4: Setting the MASCOT dynamics to Skyline.</figcaption>
</figure>

We next have to define the sampling location of the individual tips. Initially, the column **Location** should be NOT_SET for every sequence. After clicking the _Guess_ button, you can split the sequence on the vertical bar "|" again by selecting "split on character" and entering "|" in the box. However, the locations are in the fourth group, so this time choose "4" from the drop-down menu. After clicking the _OK_ button, the window should look like the one shown in the figure below:

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/TipLocations.png" alt="">
	<figcaption>Figure 5: Configuring sample locations.</figcaption>
</figure>

When leaving the priors tab and then returning to it, there will be an option to choose the population dynamics of each state separately.
Leaving the priors tab is currently necessary to make the option appear (as it forces the tab to reload), but will hopefully not be necessary in the future.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/SetSkyline.png" alt="">
	<figcaption>Figure 5: Set the Ne dynamics in both locations to Skyline dynamics.</figcaption>
</figure>


We next set the dynamics of both locations `Brazil_Northeast` and `Caribbean` to Skyline dynamics.
To the right of the skyline dynamics, we can set the number of Ne's to be estimated to 5.
This means that for each location, we will estimate 5 different effective population sizes that are equally spaced in time between the most recent sample and the root of the tree. The effective population size is estimated at five points in time. Between those five points, the effective population size is assumed to change through exponential growth.


<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/scheme.png" alt="">
	<figcaption>Figure 5: Description of the parameterization of the skyline model. In this case, 4 Ne's are estimated. Between the points where the Ne's are estimated, MASCOT-Skyline assumes exponential growth.</figcaption>
</figure>



Now, we need to set the priors for the various parameters of the model. We do this by switching to the "Priors" tab.

First, consider the effective population size parameter.  Since we have only a few samples per location, meaning little information about the different effective population sizes, we will need an informative prior. In this case we will use a log normal prior with parameters M=0 and S=1.  (These are respectively the mean and variance of the corresponding normal distribution in log space.) To use this prior, choose "Log Normal" from the dropdown menu to the right of the Ne.t:H3N2 parameter label, then click the arrow to the left of the same label and fill in the parameter values appropriately (i.e. M=0 and S=1). Ensure that the "mean in real space" checkbox remains unchecked.

The existing exponential distribution as a prior on the migration rate puts much weight on lower values while not prohibiting larger ones. For migration rates, a prior that prohibits too large values while not greatly distinguishing between very small and very *very* small values is generally a good choice. Be aware however that the exponential distribution is quite an informative prior: one should be careful that to choose a mean so that feasible rates are at least within the 95% HPD interval of the prior.  (This can be determined by clicking the arrow to the left of the parameter name and looking at the values below the graph that appears on the right.)

Finally, set the prior for the clock rate. We have a good idea about the clock rate of Influenza A/H3N2 Hemagglutinin. From previous work by other people, we know that the clock rate will be around 0.005 substitution per site per year. To include that prior knowledger, we can set the prior on the clock rate to a log normal distribution with mean in **real space**. To specify the mean in real space, make sure that the box *Mean In Real Space* is checked. If we set the S value to 0.25, we say that we expect the clock rate to be with 95% certainty between 0.00321 and 0.00731.
<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/Priors.png" alt="">
	<figcaption>Figure 6: Set up of the prior distributions.</figcaption>
</figure>


### Specify the MCMC chain length (MCMC)

Now switch to the "MCMC" tab. Here we can set the length of the MCMC
chain and decide how frequently the parameter and trees are
logged. For this dataset, 2 million iterations should be
sufficient. In order to have enough samples but not create too large
files, we can set the logEvery to 5000, so we have 401 samples
overall. Next, we have to save the `*.xml` file using _File >> Save
as_.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/MCMC.png" alt="">
	<figcaption>Figure 7: save the \*.xml.</figcaption>
</figure>

### Run the Analysis using BEAST2

Run the `*.xml` using BEAST2 or use finished runs from the *precooked-runs* folder. The analysis should take about 6 to 7 minutes. If you want to learn some more about what the migration rates we actually estimate, have a look at this blog post of Peter Beerli [http://popgen.sc.fsu.edu/Migrate/Blog/Entries/2013/3/22_forward-backward_migration_rates.html](http://popgen.sc.fsu.edu/Migrate/Blog/Entries/2013/3/22_forward-backward_migration_rates.html).

### Analyse the log file using Tracer

First, we can open the `*.log` file in tracer to check if the MCMC has converged. The ESS value should be above 200 for almost all values and especially for the posterior estimates.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LogPosterior.png" alt="">
	<figcaption>Figure 8: Check if the posterior converged.</figcaption>
</figure>

We can have a look at the marginal posterior distributions for the effective population sizes. New York is inferred to have the largest effective population size before Hong Kong and New Zealand. This tells us that two lineages that are in New Zealand are expected to coalesce quicker than two lineages in Hong Kong or New York.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LogNe.png" alt="">
	<figcaption>Figure 9: Compare the different inferred effective population sizes.</figcaption>
</figure>

In this example, we have relatively little information about the effective population sizes of each location. This can lead to estimates that are greatly informed by the prior. Additionally, there can be great differences between median and mean estimates. The median estimates are generally more reliable since they are less influence by extreme values.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/MeanMedian.png" alt="">
	<figcaption>Figure 10: Differences between mean and median estimates.</figcaption>
</figure>

We can then look at the inferred migration rates. The migration rates have the label b_migration.\*, meaning that they are backwards in time migration rates. The highest rates are from New York to Hong Kong. Because they are backwards in time migration rates, this means that lineages from New York are inferred to be likely from Hong Kong if we're going backwards in time. In the inferred phylogenies, we should therefore make the observation that lineages ancestral to samples from New York are inferred to be from Hong Kong backwards.

A more in depth explanation of what backwards migration really are can be found here [http://popgen.sc.fsu.edu/Migrate/Blog/Entries/2013/3/22_forward-backward_migration_rates.html](http://popgen.sc.fsu.edu/Migrate/Blog/Entries/2013/3/22_forward-backward_migration_rates.html)

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LogMigration.png" alt="">
	<figcaption>Figure 11: Compare the inferred migration rates.</figcaption>
</figure>

### Make the MCC tree using TreeAnnotator

Next, we want to summarize the trees. This we can do using TreeAnnotator. Open the program and then set the options as below. You have to specify the _Burnin percentage_, the _Node heights_, _Input Tree File_ and the _Output File_. After clicking _Run_ the program should summarize the trees.

<figure>
	<a id="fig:example1"></a>
	<img style="width:50%;" src="figures/TreeAnnotator.png" alt="">
	<figcaption>Figure 12: Make the maximum clade credibility tree.</figcaption>
</figure>

### Check the MCC tree using FigTree
In each logging step of the tree during the MCMC, MASCOT logs several different things. It logs the inferred probability of each node being in any possible location. In this example, these would be the inferred probabilities of being in Hong Kong, New York and New Zealand. Additonally, it logs the most likely location of each node.

After opening the MCC tree in FigTree, we can visualize several things.
To color branches, you can go to _Appearance >> Colour by_ and select *max*. This is the location that was inferred to be most often the most likely location of the node.

<figure>
<a id="fig:example1"></a>
<img style="width:100%;" src="figures/ColorsTree.png" alt="">
<figcaption>Figure 13: Inferred node locations.</figcaption>
</figure>

We can now determine if lineages ancestral to samples from New York are actually inferred to be from Hong Kong, or the probability of the root being in any of the locations.

To get the actual inferred probabilities of each node being in any of the 3 locations, you can go to _Node Labels >> Display_ an then choose Hong\_Kong, New\_York or New\_Zealand. These are the actual inferred probabilities of the nodes being in any location.

It should however be mentioned that the inference of nodes being in a particular location makes some simplifying assumptions, such as that there are no other locations (i.e. apart from the sampled locations) where lineages could have been.

Another important thing to know is that currently, we assume rates to be constant. This means that we assume that the population size of the different locations does not change over time. We also make the same assumption about the migration rates through time.

### Errors that can occur (Work in progress)

One of the errors message that can occur regularly is the following:
`too many iterations, return negative infinity`
This occurs when the integration step size of the ODE's to compute the probability of observing a phylogenetic tree in MASCOT is becoming too small.
This generally occurs if at least one migration rate is really large or at least one effective population size is really small (i.e. the coalescent rate is really high).
This causes integration steps to be extremely small, which in turn would require a lot of time to compute the probability of a phylogenetic tree under MASCOT.
Instead of doing that, this state is rejected by assigning its log probability the value negative infinity.

This error can have different origins and a likely incomplete list is the following:
1. The priors on migration rates put too much weight on really high rates. To fix this, reconsider your priors on the migration rates. Particularly, check if the prior on the migration rates make sense in comparison to the height of the tree. If, for example, the tree has a height of 1000 years, but the prior on the migration rate is exponential with mean 1, then the prior assumption is that between any two states, we expected approximately 1000 migration events.
2. The prior on the effective population sizes is too low, meaning that the prior on the coalescent rates (1 over the effective population size) is too high. This can for example occur when the prior on the effective population size was chosen to be 1/X. To fix, reconsider your prior on the effective population size.
3. There is substantial changes of the effective population sizes and/or migration rates over time that are not modeled. In that case, changes in the effective population sizes or migration rates have to be explained by population structure, which can again lead to some effective population sizes being very low and some migration rates being very high. In that case, there is unfortunately not much that can be done, since MASCOT is not an appropriate model for the dataset.
4. There is strong subpopulation structure within the different subpopulations used. In that case, reconsider if the individual sub-populations used are reasonable.


----

# Useful Links

If you interested in the derivations of the marginal approximation of the structured coalescent, you can find them here {% cite Mueller2017 --file Mascot-Tutorial/master-refs %}. This paper also explains the mathematical differences to other methods such as the theory underlying BASTA. To get a better idea of how the states of internal nodes are calculated, have a look in this paper {% cite mueller2017mascot --file Mascot-Tutorial/master-refs %}.

- MASCOT source code: [https://github.com/nicfel/Mascot](https://github.com/nicfel/Mascot)
- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Mascot-Tutorial/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users)

----

# Relevant References

{% bibliography --cited --file Mascot-Tutorial/master-refs %}
