
# make sure the working directory contains all output files from the replicate fastStructure run,
# meaning, all the .meanQ and .log files, as well as the original input file
# (use runRepStructure.pl to create a command file for replicate runs)

source("fastStructurePlotting_functions.R")
struc=read.table("nohead_6cols.str") # input structure file with population designations in second column

# collating likelihoods:
system("grep \"Marginal Likelihood =\" *.log | perl -pe 's/:Marginal Likelihood = /\t/' >likes")
likes=read.table("likes")
likes[,1]=sub(".log","",likes[,1])

means=averageBest(likelihoods=likes,top=25)
ggplotStructure(struc,means)+theme(legend.position="none")
