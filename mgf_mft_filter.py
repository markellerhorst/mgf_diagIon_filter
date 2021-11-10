# python script for the filtering of ms2 spectra in mgf format regarding potential mycofactocin species
# spectra are extracted from raw files using MSConvert with msLevel 2-2 and output format .mgf
# all spectra are scanned for the presence of fragments with 5 ppm tolerance:
# 107.0495, 136.0755, 234.1125, 235.1441, 236.1281
# written by Ellerhorst, Mark 2021
# tested on python 3.9.5

import sys
import re
import numpy as np
import pandas as pd

# check spectrum for presence of expected candidate fragments
def frag_check(entry):
    
    # calculate 5 ppm tolerance range for target fragments
    tol = 5e-06
    fragments = [107.0495, 136.0755, 236.1281]
    fragments_df = pd.DataFrame(fragments, columns = ['exact_mass'])
    fragments_df['delta_mass'] = fragments_df['exact_mass'] * tol
    fragments_df['min_mass'] = fragments_df['exact_mass'] - fragments_df['delta_mass']
    fragments_df['max_mass'] = fragments_df['exact_mass'] + fragments_df['delta_mass']
    
    # split enries by line and remove those with chars in
    entry = entry.splitlines()
    spectrum = []
    
    for line in entry:
        # check if line starts with character to remove non-spectral information
        # add each spectrum line to spectrum object
        if not bool(re.match('[a-zA-Z]', line)) and not line == '':
            spectrum.append(line)
    
    # bring spectrum into dataframe and convert str to float
    spectrum_df = pd.DataFrame([line.split(' ') for line in spectrum], columns = ['mz', 'intensity'])
    spectrum_df = spectrum_df.astype(float)
        
    # check if fragments in spectrum with 5 ppm range
    count = 0
    for index, frag in fragments_df.iterrows():
        for mz in spectrum_df.mz.values:
            if (frag['min_mass'] <= mz and frag['max_mass'] >= mz):
                count += 1
            
    # return true if at least one aglycon mass found
    if count > 0:
        return True
    
    return False
    
# add candidate precursor mass and retention time to candidate dataframe
def add_candidate(candidates, entry, cand_count):
    entry = entry.splitlines()
    candidate = {}
   
    # extract candidate mass and retention time
    for line in entry:
        if line.startswith('PEPMASS'):
            candidate['precursor_mz'] = float(line.split('=')[1].split(' ')[0])
        if line.startswith('RTINSECONDS'):
            candidate['rt_min'] = float(line.split('=')[1])/60
    # add candidate name
    candidate['candidate'] = str(cand_count)
   
    # add candidate to candidates dataframe
    candidates = candidates.append(candidate, ignore_index=True)
    
    return(candidates)
    
def main():

    # first test with one file specified as argument
    if len(sys.argv) < 2:
        sys.exit("Specify mgf file!")
    else:
        mgf_path = sys.argv[1]

    # read in mgf file and split on "END IONS"
    with open(mgf_path, 'r') as mgf_con:
        mgf_content = mgf_con.read()

        if mgf_content == '':
            sys.exit("Empty mgf provided.")
        elif not mgf_content.splitlines()[0].startswith("BEGIN IONS"):
            sys.exit("Formatting not supported.")
        else:
            entries = mgf_content.split("END IONS")
    
    # add "END IONS" to each entry again, not necessary at the moment
    for n, entry in enumerate(entries):
        entries[n] = entry + "END IONS"
      
    # scan each entry for presence of fragments:
    candidates = pd.DataFrame(columns=['precursor_mz', 'rt_min', 'candidate'])
    
    # search each spectrum in mgf for fragments of interest
    cand_count = 0
    for index, entry in enumerate(entries):
        if frag_check(entry):
            candidates = add_candidate(candidates, entries[index], cand_count)
            cand_count += 1
        
    # save candidates list if not empty
    if len(candidates) > 1:
        csv_name = str(mgf_path) + "_candidates.csv"
        candidates.to_csv(csv_name, index=False)
    else:
        print("No candidates found.")
    
if __name__ == '__main__':
    main()
