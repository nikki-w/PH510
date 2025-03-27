import MonteCarlo_class as MC

MC_rand = MC.MC(D=3, N=100000, seed=134)

points = MC_rand.points()

print(points)
