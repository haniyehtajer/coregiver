import os

def generate_input_uniform(num_particles, particle_mass, particle_cmf, output_filename):
    """
    Generates a simplified DBCT input file with particle hash, mass, and CMF.

    Args:
        num_particles (int): The total number of particles to generate.
        particle_mass (float): The mass assigned to each particle.
        particle_cmf (float): The core mass fraction (CMF) assigned to each particle.
        output_filename (str): The name of the output file.
    """
    print(f"Generating {num_particles} particles with mass={particle_mass} and CMF={particle_cmf}...")
    print(f"Output file: {output_filename}")

    try:
        # Open the file in write mode. 'w' will create the file if it doesn't exist,
        # or overwrite it if it does. Use 'a' if you want to append to an existing file.
        with open(output_filename, 'w') as f:
            # Loop to generate and write particle data
            for i in range(num_particles):
                # Hash generation: Similar to the C code, starting from 1
                # For simplicity, we're treating all as one type of particle for hashing.
                particle_hash = i + 1

                # Write data in tab-separated format
                # Using f-strings for easy formatting and precision for floats
                f.write(f"{particle_hash}\t{particle_mass:.8e}\t{particle_cmf:.8e}\n")

        print(f"Successfully generated '{output_filename}' with {num_particles} entries.")

    except IOError as e:
        print(f"Error writing to file '{output_filename}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def generate_input_step_function(num_particles, particle_mass, semi_majors, inner_cmf, outer_cmf, boundary, output_filename):
    """
    Generates a simplified DBCT input file with particle hash, mass, and CMF.

    Args:
        num_particles (int): The total number of particles to generate.
        particle_mass (float): The mass assigned to each particle.
        semi_majors (array): Semi-major axes of the particles.
        inner_cmf (float): Desired CMF for the particles inside the CMF change boundary
        outer_cmf (float): Desired CMF for the particles outside the CMF change boundary
        boundary (float): The CMF change boundary
        output_filename (str): The name of the output file.
    """
    print(f"Output file: {output_filename}")

    try:
        # Open the file in write mode. 'w' will create the file if it doesn't exist,
        # or overwrite it if it does. Use 'a' if you want to append to an existing file.

        with open(output_filename, 'w') as f:
            # Loop to generate and write particle data
            for i in range(num_particles):
                # Hash generation, starting from 1
                particle_hash = i + 1
                semi = semi_majors[i]
                if semi <= boundary:
                    particle_cmf = inner_cmf
                else:
                    particle_cmf = outer_cmf

                # Write data in tab-separated format
                # Using f-strings for easy formatting and precision for floats

                f.write(f"{particle_hash}\t{particle_mass:.8e}\t{particle_cmf:.8e}\n")

        print(f"Successfully generated '{output_filename}' with {num_particles} entries.")

    except IOError as e:
        print(f"Error writing to file '{output_filename}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def get_inner_outer_mass(m_emb, m_pl, boundary, a_emb, a_pl):
    m_inner = 0
    m_outer = 0
    for a in a_pl:
        if a <= boundary:
            m_inner = m_inner + m_pl
        else:
            m_outer = m_outer + m_pl
    
    for a in a_emb:
        if a <= boundary:
            m_inner = m_inner + m_emb
        else:
            m_outer = m_outer + m_emb
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