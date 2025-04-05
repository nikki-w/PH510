"""MIT License

Copyright (c) [2025] [Nikki Walker]

File that uses MC class"""
import montecarlo_class as MC


# Test Case 1: n-dimension sphere

# Initialise MC in 2D
MC_rand_2D = MC.MC(d=2, n=100000000, seed=134)

# Calculate exact volume for D case
volume_2D_exact = MC_rand_2D.volume_exact()

# Calculate estimated volume using MC
volume_2D_MC = MC_rand_2D.volume(r=1)
volume_2D_MC_mean, volume_2D_MC_err = MC_rand_2D.statistics(volume_2D_MC)

if MC_rand_2D.rank == 0:
    print("Test Case 1: n-dimension sphere", "\n")
    print("2D sphere:", "\n")
    print(f"Exact volume: {volume_2D_exact:.6f}")
    print(f"MC estimate: {volume_2D_MC:.6f} +/- {volume_2D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_2D_MC - volume_2D_exact)/volume_2D_exact*100:.2f}% \n",
 "\n")

# Initialise MC in 3D
MC_rand_3D = MC.MC(d=3, n=100000000, seed=134)

# Calculate exact volume for D case
volume_3D_exact = MC_rand_3D.volume_exact()

# Calculate estimated volume using MC
volume_3D_MC = MC_rand_3D.volume(r=1)
volume_3D_MC_mean, volume_3D_MC_err = MC_rand_3D.statistics(volume_3D_MC)

if MC_rand_3D.rank == 0:
    print("3D sphere:", "\n")
    print(f"Exact volume: {volume_3D_exact:.6f}")
    print(f"MC estimate: {volume_3D_MC:.6f} +/- {volume_3D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_3D_MC - volume_3D_exact)/volume_3D_exact*100:.2f}% \n",
 "\n")

# Initialise MC in 4D
MC_rand_4D = MC.MC(d=4, n=100000000, seed=134)

# Calculate exact volume for D case
volume_4D_exact = MC_rand_4D.volume_exact()

# Calculate estimated volume using MC
volume_4D_MC = MC_rand_4D.volume(r=1)
volume_4D_MC_mean, volume_4D_MC_err = MC_rand_4D.statistics(volume_4D_MC)

if MC_rand_4D.rank == 0:
    print("4D sphere:", "\n")
    print(f"Exact volume: {volume_4D_exact:.6f}")
    print(f"MC estimate: {volume_4D_MC:.6f} +/- {volume_4D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_4D_MC - volume_4D_exact)/volume_4D_exact*100:.2f}% \n",
 "\n")

# Initialise MC in 5D
MC_rand_5D = MC.MC(d=5, n=100000000, seed=134)

# Calculate exact volume for D case
volume_5D_exact = MC_rand_5D.volume_exact()

# Calculate estimated volume using MC
volume_5D_MC = MC_rand_5D.volume(r=1)
volume_5D_MC_mean, volume_5D_MC_err = MC_rand_5D.statistics(volume_5D_MC)

if MC_rand_5D.rank == 0:
    print("5D sphere:", "\n")
    print(f"Exact volume: {volume_5D_exact:.6f}")
    print(f"MC estimate: {volume_5D_MC:.6f} +/- {volume_5D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_5D_MC - volume_5D_exact)/volume_5D_exact*100:.2f}% \n",
 "\n")
