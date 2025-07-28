import os
import numpy as np

def generate_input_uniform(num_particles, particle_mass, particle_cmf, output_filename):
    """
    Generates a simplified DBCT input file with particle hash, mass, and CMF.

    Args:
        num_particles (int): The total number of particles to generate.
        particle_mass (float or array-like): Mass(es) for the particles. Can be a single float or an array of floats.
        particle_cmf (float): The core mass fraction (CMF) assigned to each particle.
        output_filename (str): The name of the output file.
    """
    particle_mass = np.array(particle_mass)

    if particle_mass.size != 1 and particle_mass.size != num_particles:
        raise ValueError(f"particle_mass must be a single value or an array of length {num_particles}")

    print(f"Generating {num_particles} particles with CMF={particle_cmf}...")
    print(f"Output file: {output_filename}")

    try:
        with open(output_filename, 'w') as f:
            for i in range(num_particles):
                particle_hash = i + 1
                mass = particle_mass if particle_mass.size == 1 else particle_mass[i]
                f.write(f"{particle_hash}\t{mass:.8e}\t{particle_cmf:.8e}\n")

        print(f"Successfully generated '{output_filename}' with {num_particles} entries.")

    except IOError as e:
        print(f"Error writing to file '{output_filename}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")



def generate_input_step_function(num_particles, particle_mass, semi_majors,
                                 inner_cmf, outer_cmf, boundary, output_filename):
    """
    Generates a simplified DBCT input file with particle hash, mass, and CMF based on a step-function CMF.

    Args:
        num_particles (int): The total number of particles to generate.
        particle_mass (float or array-like): The mass(es) assigned to each particle.
        semi_majors (array-like): Semi-major axes of the particles.
        inner_cmf (float): Desired CMF for particles inside the CMF change boundary.
        outer_cmf (float): Desired CMF for particles outside the CMF change boundary.
        boundary (float): The CMF change boundary.
        output_filename (str): The name of the output file.
    """
    particle_mass = np.array(particle_mass)
    semi_majors = np.array(semi_majors)

    if semi_majors.size != num_particles:
        raise ValueError("semi_majors must have the same length as num_particles.")

    if particle_mass.size != 1 and particle_mass.size != num_particles:
        raise ValueError("particle_mass must be a scalar or an array with length equal to num_particles.")

    print(f"Output file: {output_filename}")

    try:
        with open(output_filename, 'w') as f:
            for i in range(num_particles):
                particle_hash = i + 1
                semi = semi_majors[i]
                mass = particle_mass if particle_mass.size == 1 else particle_mass[i]

                particle_cmf = inner_cmf if semi <= boundary else outer_cmf

                f.write(f"{particle_hash}\t{mass:.8e}\t{particle_cmf:.8e}\n")

        print(f"Successfully generated '{output_filename}' with {num_particles} entries.")

    except IOError as e:
        print(f"Error writing to file '{output_filename}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")



def get_inner_outer_mass(m_emb, m_pl, boundary, a_emb, a_pl):
    m_inner = 0
    m_outer = 0

    m_emb = np.array(m_emb)
    m_pl = np.array(m_pl)

    for mass, a in zip(m_pl, a_pl):
        if a <= boundary:
            m_inner += mass
        else:
            m_outer += mass

    for mass, a in zip(m_emb, a_emb):
        if a <= boundary:
            m_inner += mass
        else:
            m_outer += mass

    return m_inner, m_outer


def get_inner_outer_cmf(cmf_initial, x, m_inner, m_outer):
    m_total = m_inner + m_outer
    iron_moved = x * (cmf_initial * m_outer)
    cmf_inner_new = (cmf_initial * m_inner + iron_moved)/m_inner
    cmf_outer_new = (cmf_initial * m_outer - iron_moved)/m_outer

    m_Fe_init = cmf_initial * m_total
    m_Fe_final = cmf_inner_new * m_inner + cmf_outer_new * m_outer
    diff = m_Fe_final - m_Fe_init
    if diff > 1e-10:
        print("Warning: difference in total Iron Mass is more than 1e-10.")
        
    return cmf_inner_new, cmf_outer_new

def get_dual_cmf(m_emb, m_pl, semis_emb, semis_pl, boundary, x, cmf_init): #mass of emb and pl, semi major of emb and pl, boundary between
                                                                           #inner and outer disk, fraction of iron transport, initial cmf
    m_in, m_out = get_inner_outer_mass(m_emb, m_pl, boundary, semis_emb, semis_pl)
    cmf_in, cmf_out = get_inner_outer_cmf(cmf_init, x, m_in, m_out)
    return cmf_in, cmf_out