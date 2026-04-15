import numpy as np

# -----------------------------
# PARAMETERS
# -----------------------------
N_particles = 5000

mu_true = 1.0
mu_biased = 0.3

shield_thickness = 5.0

w_min = 0.01
p_survive = 0.1

# interaction probabilities
p_absorb = 0.3
p_scatter = 0.5
p_pair = 0.2

# -----------------------------
# RESULTS
# -----------------------------
E_out_total = 0.0
E_abs_total = 0.0
E_sim_total = 0.0

# -----------------------------
# MAIN LOOP
# -----------------------------
for i in range(N_particles):

    # stack of particles
    particles = []

    # initial photon
    particles.append({
        "E": 1.5,   # must allow pair production
        "x": 0.0,
        "dir": 1.0, # 1D direction
        "w": 1.0
    })

    E_sim_total += 1.5

    # -----------------------------
    # TRACK ALL PARTICLES
    # -----------------------------
    while particles:

        p = particles.pop()

        E = p["E"]
        x = p["x"]
        direction = p["dir"]
        w = p["w"]

        alive = True

        while alive:

            # -----------------------------
            # SAMPLE FREE PATH (biased)
            # -----------------------------
            r = np.random.rand()
            s = -np.log(r) / mu_biased

            d_out = (shield_thickness - x) if direction > 0 else x

            # -----------------------------
            # WEIGHT UPDATE
            # -----------------------------
            weight_factor = (mu_true / mu_biased) * np.exp(-(mu_true - mu_biased) * s)
            w *= weight_factor

            # -----------------------------
            # ESCAPE?
            # -----------------------------
            if s >= d_out:
                E_out_total += w * E
                break

            # move particle
            x += direction * s

            # -----------------------------
            # INTERACTION TYPE
            # -----------------------------
            xi = np.random.rand()

            if xi < p_absorb:
                # absorption
                E_abs_total += w * E
                break

            elif xi < p_absorb + p_scatter:
                # scattering (toy)
                E *= 0.8

            else:
                # -----------------------------
                # PAIR PRODUCTION
                # -----------------------------
                # photon disappears → create 2 new photons

                new_E = 0.511

                particles.append({
                    "E": new_E,
                    "x": x,
                    "dir": 1.0,
                    "w": w
                })

                particles.append({
                    "E": new_E,
                    "x": x,
                    "dir": -1.0,
                    "w": w
                })

                # parent photon terminated
                break

            # -----------------------------
            # RUSSIAN ROULETTE
            # -----------------------------
            if w < w_min:
                if np.random.rand() < p_survive:
                    w /= p_survive
                else:
                    break

# -----------------------------
# RESULTS
# -----------------------------
print("E_simulated =", E_sim_total)
print("E_out       =", E_out_total)
print("E_abs       =", E_abs_total)
print("E_out + E_abs =", E_out_total + E_abs_total)