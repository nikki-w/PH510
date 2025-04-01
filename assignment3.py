"""MIT License

Copyright (c) [2025] [Nikki Walker]

File that uses MC class"""
import montecarlo_class as MC

# Initialise MC in 2D
MC_rand = MC.MC(d=2, n=1000000000, seed=134)

# Calculate exact volume for D case 
volume_2D_exact = MC_rand.volume_exact()

# Calculate estimated volume using MC
volume_2D_MC = MC_rand.volume(r=1)

if MC_rand.rank == 0:
    print(f"Exact volume: {volume_2D_exact:.6f}")
    print(f"MC estimate: {volume_2D_MC:.6f}") 
