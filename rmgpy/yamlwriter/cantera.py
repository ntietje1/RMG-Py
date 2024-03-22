from datetime import datetime
from rmgpy.chemkin import get_species_identifier



def write_cantera(
    spcs,
    rxns,
    surface_site_density=None,
    solvent=None,
    solvent_data=None,
    path="chem.yml",
):
    """
    Writes yaml file depending on the type of system (gas-phase, catalysis).
    Writes beginning lines of yaml file, then uses yaml.dump(result_dict) to write species/reactions info. 
    """
    
    """
    -- PHASES:
    ---- BLOCK ONE LIKE THIS (hard coded):

phases:
- name: gas
  thermo: ideal-gas
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]
  
    ---- BLOCK TWO LIKE THIS (depends on spcs list):
    
species: [Ar, He, Ne, N2, octane(1), oxygen(2), 'CC[CH]CCCCC(9)', 'CCC[CH]CCCC(12)', 'C[CH]CCCCCC(15)', '[CH2]CCCCCCC(16)', '[O]O(17)', OO(20), C8H16(28), C8H17O2(43), CCCCCCC(C)OO(44), C8H16(53), C8H16(54), C8H17O2(67), CCCCC(CCC)OO(68), C8H17O2(84), CCCCCC(CC)OO(85), 'CCCCCCCCO[O](99)', CCCCCCCCOO(100)]
kinetics: gas
    
    ---- BLOCK THREE LIKE THIS (hard coded):
    
transport: mixture-averaged
state: {T: 300.0, P: 1 atm}

    -- ELEMENTS:
    ---- LIKE THIS (hard coded):
    
elements:
- symbol: Ci
  atomic-weight: 13.003
- symbol: D
  atomic-weight: 2.014
- symbol: Oi
  atomic-weight: 17.999
- symbol: T
  atomic-weight: 3.016
- symbol: X
  atomic-weight: 195.083
  
    --- BODY: YAML
    --- LIKE THIS:
    
    species:
- name: Ar
  composition: {Ar: 1}
  thermo:
    model: NASA7
    temperature-ranges: [100.0, 4879.8, 5000.0]
    data:
    - [2.5, -3.01680531e-12, 3.74582141e-15, -1.50856878e-18, 1.86626471e-22, -774.692768, 3.18514899]
    - [4.28461071, -1.45494649e-03, 4.44804306e-07, -6.04359642e-11, 3.07921551e-15, -2525.81822, -8.26300883]
    note: |-
      Thermo library corrected for liquid phase: primaryThermoLibrary + Solvation correction with octane as solvent and solute estimated using Solute
      library: Ar
  transport:
    model: gas
    geometry: atom
    well-depth: 136.5
    diameter: 3.33
    note: NOx2018
  note: Ar
  - <<< REPEAT FOR EACH SPECIES >>>>
  
reactions:
- equation: '[O]O(17) + C[CH]CCCCCC(15) <=> oxygen(2) + octane(1)'  # Reaction 1
  rate-constant: {A: 0.02759059, b: 3.802, Ea: 6.324}
  note: |-
    Reaction index: Chemkin #1; RMG #13
    Template reaction: H_Abstraction
    Flux pairs: [O]O(17), oxygen(2); C[CH]CCCCCC(15), octane(1);
    Estimated using template [X_H;C_rad/H/NonDeC] for rate rule [Orad_O_H;C_rad/H/NonDeC]
    Euclidian distance = 2.0
    family: H_Abstraction
- <<< REPEAT FOR EACH REACTION >>>>
  
  OUTPUT IN THIS ORDER:
  
  <PREFACE>
  <PHASES>
  <TRANSPORT>
  <ELEMENTS>
  <SPECIES>
  
  <BODY>
  
    """

def is_surface(spcs):
    """
    Returns True if any species in the list contains a surface site.
    """
    for spc in spcs:
        if spc.contains_surface_site():
            return True
    return False

def write_preface() -> str:
    """
    Writes the preface block of the yaml file.
    """
    
    lines = [] # use list for fast string concatenation
    
    # generator line
    lines.append("generator: RMG\n")

    # datetime object containing current date and time
    now = datetime.now()
    dt_string = now.strftime("%a, %d %b %Y %H:%M:%S")
    lines.append(f"date: {dt_string}\n")

    # units line
    lines.append("\nunits: {length: cm, time: s, quantity: mol, activation-energy: kcal/mol}\n")
    lines.append("\n")
    
    return ''.join(lines) # return as string (could also return as list here)

def write_elements_block():
    """
    todo:
    """
    return """
elements:
- symbol: Ci
  atomic-weight: 13.003
- symbol: D
  atomic-weight: 2.014
- symbol: Oi
  atomic-weight: 17.999
- symbol: T
  atomic-weight: 3.016
- symbol: X
  atomic-weight: 195.083
"""

def write_species_block(spcs):
    """
    todo:
    """
    sorted_species = sorted(spcs, key=lambda spc: spc.index)
    return "species: [" + ', '.join(f"'{get_species_identifier(spc)}'" for spc in sorted_species) + "]"

def write_nonsurface_species(spcs) -> str:
    """
    todo:
    """
    block1 = """
phases:
- name: gas
  thermo: ideal-gas
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]"""
  
    block2 = write_species_block(spcs) + "\nkinetics: gas"
    
    block3 = """
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}
        """
        
    block4 = write_elements_block()

    return block1, block2, block3, block4

def write_surface_species(spcs, rxns, surface_site_density) -> str:
    """
    todo:
    """
    surface_spcs = [spc for spc in spcs if spc.contains_surface_site()]
    gas_spcs = [spc for spc in spcs if not spc.contains_surface_site()]

    block1 = f"""
phases:
- name: gas
  thermo: ideal-gas
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]
  {write_species_block(gas_spcs)}
  kinetics: gas
  reactions: [gas_reactions]"""

    block2 = """
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}"""

    block3 = f""" 
- name: {surface_spcs[0].smiles.replace("[","").replace("]","")}_surface
  thermo: ideal-surface
  adjacent-phases: [gas]
  elements: [H, D, T, C, Ci, O, Oi, N, Ne, Ar, He, Si, S, F, Cl, Br, I, X]
  {write_species_block(surface_spcs)}
  kinetics: surface
  reactions: [surface_reactions]     
  site-density: {surface_site_density * 1e-4 }
        """

    block4 = write_elements_block()

    return block1, block2, block3, block4
    
    
def write_body(spcs, rxns, surface_site_density) -> str:
    """
    Writes the body of the yaml file.
    """
    
    if is_surface(spcs):
        # todo: WRITE SURFACE SPECIES
        # result_dict = get_mech_dict_surface(
        #     spcs, rxns, solvent=solvent, solvent_data=solvent_data
        # )
        block1, block2, block3, block4 = write_surface_species(
            spcs, rxns, surface_site_density
        )
    else:
        # todo: WRITE NON-SURFACE SPECIES
        # # get_mech_dict writes yaml files without creating separate
        # result_dict = get_mech_dict_nonsurface(
        #     spcs, rxns, solvent=solvent, solvent_data=solvent_data
        # )
        block1, block2, block3, block4 = write_nonsurface_species(spcs)
