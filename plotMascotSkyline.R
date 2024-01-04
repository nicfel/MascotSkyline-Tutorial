# Install packages if they are not already installed
if (!requireNamespace("ggtree", quietly = TRUE)) install.packages("ggtree")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("treeio", quietly = TRUE)) install.packages("treeio")
if (!requireNamespace("coda", quietly = TRUE)) install.packages("coda")


# Load the necessary libraries
library(ggtree)
library(ape)
library(ggplot2)
library(treeio)
library(coda)

# Define the Okabe-Ito color scale (8 colors)
okabe_ito_scale <- c("Caribbean"="#0072B2", "Brazil_Northeast"="#D55E00")

# Read in the BEAST MCC tree, after adding the path to the mcc file
mcc_tree <- read.beast("/Users/nmueller/Library/CloudStorage/OneDrive-FredHutchinsonCancerCenter/Documents/github/MascotSkyline-Tutorial/precooked_runs/ZIKV.mcc.trees")  # Change to your file path

# Plot the tree using ggtree
p_tree = ggtree(mcc_tree, aes(color=max)) +
  scale_color_manual(values = okabe_ito_scale) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
  )+
  geom_tippoint(aes(color = max))  # Map the tip points' color to the location
plot(p_tree)
# ggsave(plot=p_tree, width=4, height=4, "tree.pdf")


# plot the effective population sizes over time
log = read.table(header=T, sep="\t", 
  "/Users/nmueller/Library/CloudStorage/OneDrive-FredHutchinsonCancerCenter/Documents/github/MascotSkyline-Tutorial/precooked_runs/ZIKV.log"
  )

# take a 10% burnin
burnin = 0.1
burnin = round(burnin*length(log$Sample))
log = log[(burnin+1):length(log$Sample),]

# get the minimum tree height
min_tree_height = min(log$Tree.height)

# use 100 points to interpolate the Ne dynamics 
times_points_interpl = seq(0,min_tree_height,length.out=100)

# define the time of the most recent sampled individual
mrsi = 2016.762295

# get all the headers of the log file
param_labels = labels(log)[[2]]

# get the locations by checking all the names in SkykineNe.locationName.
all_locations = c()
for (i in seq(1, length(param_labels))){
  if (startsWith(param_labels[i], "SkylineNe.")){
    all_locations = c(all_locations, strsplit(param_labels[i], ".", fixed=T)[[1]][2])
  }
}
locations = unique(all_locations)

# initialize data fram
dat = data.frame()

for (i in seq(1, length(locations))){
  # get the parameters for the location
  params = param_labels[startsWith(param_labels,
                      paste("SkylineNe.", locations[[i]], ".",sep=""))]
  tmp.ne = c()
  # get the time points where the Ne's where estimated by getting length(params) points
  # between 0 and 1
  change_points = seq(0,1,length.out=length(params))
  # for each iteration in the posterior after the burnin, get the Ne's at the time
  # points used for plotting
  for (j in seq(1, length(log$Sample))){
    times = change_points*log[j, "Tree.height"]
    vals = log[j, params]
    tmp.ne = c(tmp.ne, approx(times, vals, xout=times_points_interpl, method = "linear")$y)
  }
  # Now, we have a grid of Ne points and can interpolate the Ne's and compute
  # 95% HPD intervals
  for (j in seq(1, length(times_points_interpl))){
    # get the exponential of all values (the Ne's are logged in log space)
    vals=exp(as.mcmc(tmp.ne[seq(j, length(tmp.ne), length(times_points_interpl))]))
    hpd.5 = HPDinterval(vals, prob=0.5)
    hpd.95 = HPDinterval(vals, prob=0.95)
    timestart = mrsi-times_points_interpl[[j]]
    dat = rbind(dat, data.frame(time = timestart,
                                            l.5=hpd.5[1,"lower"], u.5=hpd.5[1,"upper"],
                                            l.95=hpd.95[1,"lower"], u.95=hpd.95[1,"upper"],
                                            location=locations[[i]],
                                            method="skygrid"))
  }
}

# plot the skylines
p_skyline = ggplot(data=dat, aes(x=time, fill=location))+
      geom_ribbon(aes(ymin=l.5, ymax=u.5), alpha=0.5)+ # plots the 50% HPD
      geom_ribbon(aes(ymin=l.95, ymax=u.95), alpha=0.5)+ # plots the 95% HPD
      scale_fill_manual(values=okabe_ito_scale)

plot(p_skyline)
      
