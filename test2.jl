using StructuredLight #Loads the package

rs = LinRange(-5,5,256) #Define a linear grid of points
zs = LinRange(0,1,32)

E0 = lg(rs,rs,zs) #Calculates the fundamental Laguerre-Gaussian mode

show_animation(E0) #visualizes the mode

E = free_propagation(E0,rs,rs,1) #Propagates the mode through a distance of z=1

visualize(E) #visualizes the evolved mode