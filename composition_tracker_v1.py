# composition_tracker.py

import numpy as np
import sys
from astropy import units as u
import pandas as pd

# User defines (can be set as default arguments or passed in)
user_core_density = 7874.0  # kg/m^3
user_mantle_density = 3000.0  # kg/m^3

# Units
kg_per_m_cubed = u.kg / u.m**3
Msun_per_au_cubed = u.Msun / u.au**3

# Densities, unit conversion from kg/m^3 to M_sun/AU^3 to match simulation units
core_density = (user_core_density * kg_per_m_cubed).to(Msun_per_au_cubed)
mantle_density = (user_mantle_density * kg_per_m_cubed).to(Msun_per_au_cubed)


def organize_compositions(init_compositions):
    """Data organizing function

    Takes in raw string data extracted from composition input file
    Creates a list filled with properties of each particle
    The properties include hash, mass, and CMF

    Parameters:
    init_compositions (list) -- raw particle data from input file

    Returns:
    compositions (list) -- nested list with properly formatted particle data
    """
    compositions = []
    for particle in init_compositions:
        particle_data = []
        for i in range(len(particle)):
            if i == 0:
                particle_data.append(int(particle[i]))  # The particle's REBOUND hash
            elif i < 3:
                particle_data.append(float(particle[i]))  # Adds its mass than CMF
        compositions.append(particle_data)
    return compositions


def calc_radius(m, cmf, mantle_density, core_density):
    """Return the radius of an object based on its CMF and mass"""
    m_core = cmf * m
    m_mantle = m - m_core
    rho_tot = (m_core + m_mantle) / ((m_core / core_density) + (m_mantle / mantle_density))
    v_tot = m / rho_tot
    r_tot = (3 / (4 * np.pi) * v_tot) ** (1 / 3)
    return r_tot


def track_composition(collision_report_file, composition_input_file, ejection_file):
    """Main function

    Tracks how collisions change the CMFs of objects from REBOUND fragmentation sim
    Returns the final compositions of all remaining objects
    Has 2 main sections for mergers and disruptive collisions (accretive and erosive)


    Parameters:
    collision_report_file (str) -- pathway to collision report file
    composition_input_file (str) -- pathway to the DBCT input file
    ejection_file (str) -- pathway to file that lists objects ejected from sim

    Returns:
    compositions (list) -- nested list with compositional data of final objects
    """
    # Extracts compositional data
    f = open(composition_input_file, 'r')
    raw_compositions = [line.split() for line in f.readlines()]
    compositions = organize_compositions(raw_compositions)  # This list keeps track of all the particle data
    for obj in compositions:
        if obj[2] > 1.0 or obj[2] < 0.0:
            print('ERROR: CMF does not have a realistic value')
            sys.exit(1)
    f.close()

    # Extracts collisional data
    f = open(collision_report_file, 'r')
    collision_blocks = f.read().split("\n")
    collisions = [block for block in collision_blocks if len(block) > 0]
    f.close()

    destroyed_object_hashes = []  # Keeps track of objects that get destroyed in collisions

    ############### START OF MAIN LOOP ###############
    # Iterates through every collision from sim
    # Determines the change in objects' compositions from each collision
    # Adds new fragments to the compositions list when necessary
    for i in range(len(collisions)):
        collision = collisions[i].split()
        time = float(collision[0])
        collision_type = int(collision[1])
        if collision_type == 0:  # Elastic bounces do nothing so loop moves to next collision
            continue
        sim_impact_param = float(collision[2])
        target_hash = int(collision[3])
        largest_remnant_mass = float(collision[4])  # Mass of the target after collision
        target_sim_radius = float(collision[5]) / 5
        proj_hash = int(collision[6])
        proj_sim_radius = float(collision[7]) / 5
        no_frags = int((len(collision) - 8) / 2)
        frag_hashes = [int(collision[j * 2 + 6]) for j in range(1, no_frags + 1)]
        frag_masses = [float(collision[j * 2 + 7]) for j in range(1, no_frags + 1)]

        # Determines if the projectile collided with something not in the compositions list (Star, Jupiter, etc.)
        big_obj_collision_flag = 1
        for obj in compositions:
            if target_hash == obj[0]:
                big_obj_collision_flag += -1
                break
        if big_obj_collision_flag == 1:
            destroyed_object_hashes.append(proj_hash)  # Projectile is always considered destroyed in this type of collision
            if no_frags != 0:  # If fragments created, adds those to compositions list with CMF of 0.3
                for j in range(no_frags):
                    frag_data = [frag_hashes[j], frag_masses[j], 0.3]
                    compositions.append(frag_data)
            continue  # Moves to the next collision

        targ_idx = [idx for idx in range(len(compositions)) if int(compositions[idx][0]) == target_hash][0]  # Index of target in the compositions list
        proj_idx = [idx for idx in range(len(compositions)) if int(compositions[idx][0]) == proj_hash][0]  # Index of projectile in the compositions list
        target_mass = float(compositions[targ_idx][1])
        proj_mass = float(compositions[proj_idx][1])
        target_core_frac = compositions[targ_idx][2]
        proj_core_frac = compositions[proj_idx][2]

        # Total core mass available
        total_core_mass = target_core_frac * target_mass + proj_core_frac * proj_mass

        # Total mantle mass available
        total_mantle_mass = target_mass + proj_mass - total_core_mass

        # This sequence estimates the outer radii of the target and projectile
        target_radius = calc_radius(target_mass, target_core_frac, mantle_density, core_density)
        proj_radius = calc_radius(proj_mass, proj_core_frac, mantle_density, core_density)
        sine_of_impact_angle = sim_impact_param / (target_sim_radius + proj_sim_radius)
        impact_parameter = (target_radius + proj_radius) * sine_of_impact_angle  # New impact parameter with updated radii

    ############### PERFECT MERGER ###############
        # basically weighted average of initial target compoisition and mass with the projectile composition and mass
        if collision_type == 1:  # perfect merger
            compositions[targ_idx][2] = ((target_core_frac * target_mass) + (proj_core_frac * proj_mass)) / largest_remnant_mass

    ############### DISRUPTIVE COLLISIONS ###############
        elif collision_type == 2 or collision_type == 3 or collision_type == 4:
            # Collision type 2 -> partial accretion
            if collision_type == 2:
                # target has accreted mass, if its mass is larger than total core mass available,
                # Then all core mass will be target's core, and fragments will be mantle material only
                if largest_remnant_mass > total_core_mass:
                    CMF_lr = total_core_mass / largest_remnant_mass
                    CMF_frag = 0

                # If Mlr is less than total core mass available, then Mlr will be pure Iron,
                # and fragments will share the rest of core material equally
                else:
                    CMF_lr = 1
                    CMF_frag = (total_core_mass - largest_remnant_mass) / (target_mass + proj_mass - largest_remnant_mass)

                compositions[targ_idx][2] = CMF_lr  # Changes target CMF to largest remnant CMF for accretion

            # Collision type 3 and 4 -> erosive collisions
            else:
                # If Mlr is larger than core mass available, then all core is assigned to Mlr
                # All frags will be mantle material only
                if largest_remnant_mass > total_core_mass:
                    CMF_lr = total_core_mass / largest_remnant_mass
                    CMF_frag = 0

                # If Mlr is less than core mass available, then its CMF will be 1
                # and fragments will share core material equally
                else:
                    CMF_lr = 1
                    CMF_frag = (total_core_mass - largest_remnant_mass) / (target_mass + proj_mass - largest_remnant_mass)

                compositions[targ_idx][2] = CMF_lr  # Changes target CMF to largest remnant CMF for erosion

            # Target CMF calculation tolerance
            if compositions[targ_idx][2] - 1.0 > 0.0 and compositions[targ_idx][2] - 1.0 < 1.0e-5:  # If target CMF is just above 1.0
                compositions[targ_idx][2] = 1.0
            if compositions[targ_idx][2] < 0.0 and compositions[targ_idx][2] > -1.0e-5:  # If target CMF is just below 0.0
                compositions[targ_idx][2] = 0.0
          
            frag_core_frac = CMF_frag  # All fragments get the same CMF

            # Fragment CMF calculation tolerance
            if frag_core_frac - 1.0 > 0.0 and frag_core_frac - 1.0 < 1.0e-5:  # If fragment CMF is just above 1.0
                frag_core_frac += 1.0 - frag_core_frac
            if frag_core_frac < 0.0 and frag_core_frac > -1.0e-5:  # If fragment CMF is just below 0.0
                frag_core_frac += 0.0 - frag_core_frac

            # Fragment CMF error
            if frag_core_frac < 0.0:
                print('ERROR: Fragment CMF is negative at time: ', time)
                sys.exit(1)
            elif frag_core_frac > 1.0:
                print('ERROR: Fragment CMF is greater than 1.0 at time and is : ', time, frag_core_frac)
                sys.exit(1)

            for j in range(no_frags):
                if frag_masses[j] < 0:
                    print('ERROR: Fragment mass is negative at time: ', time)
                    sys.exit(1)
                else:
                    frag_data = [frag_hashes[j], frag_masses[j], frag_core_frac]
                    compositions.append(frag_data)  # Fragment now added to the object tracking list

        compositions[targ_idx][1] = largest_remnant_mass  # Changes target mass to largest remnant mass

        # Checks to see if there are any errors with the largest remnant's properties
        for j in range(len(compositions[targ_idx])):
            if j == 1:
                if compositions[targ_idx][j] < 0.0:
                    print('ERROR: Largest remnant mass is negative at time: ', time)
                    sys.exit(1)
            elif j > 1:
                if compositions[targ_idx][j] < 0.0:
                    print('ERROR: Largest remnant CMF is negative at time: ', time)
                    sys.exit(1)
                elif compositions[targ_idx][j] > 1.0:
                    print('ERROR: Largest remnant CMF is greater than 1.0 at time: ', time)
                    sys.exit(1)

        # Checks to make sure the projectile isn't a second largest remnant before deletion
        for hsh in frag_hashes:
            if hsh != proj_hash and hsh == frag_hashes[-1]:
                destroyed_object_hashes.append(proj_hash)
            elif hsh == proj_hash:
                break
            else:
                continue

        if collision_type == 1:
            destroyed_object_hashes.append(proj_hash)  # The projectile in a merger is always destroyed

    ############### END OF MAIN LOOP ###############

    # Removes destroyed objects from compositions list
    # Need to iterate backwards when popping to avoid index issues
    for hsh in destroyed_object_hashes:
        for i in range(len(compositions) - 1, -1, -1):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break

    # Removes ejected objects from compositions list
    f = open(ejection_file, 'r')
    ejections_raw = [line.split() for line in f.readlines()]
    ejections = [abs(int(ejections_raw[i][1])) for i in range(len(ejections_raw))]
    f.close()

    # Need to iterate backwards when popping to avoid index issues
    for hsh in ejections:
        for i in range(len(compositions) - 1, -1, -1):
            if compositions[i][0] == hsh:
                compositions.pop(i)
                break
    return compositions


def write_output(compositions, composition_output_file):
    """Writes final objects and their propeties to output file"""
    with open(composition_output_file, "w") as f:
        for obj in compositions:
            for item in obj:
                f.write(str(item) + ' ')
            f.write('\n')