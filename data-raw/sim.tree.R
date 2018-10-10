# Save the simulated tree in an object of class phylo

# Read simulated tree with ape
sim.tree <- ape::read.tree(text="((((A:0.1,B:0.1):0.2,(F:0.1,C:0.2):0.1):0.5,H:0.1):0.2,(D:0.7,(G:0.2,E:0.5):0.2):0.3);")

# Save simulated tree in rda format
devtools::use_data( sim.tree )
