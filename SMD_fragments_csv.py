import itertools
from pyopenms import *
import pandas as pd

def SMD_fragments(list_of_peptides, lower_mass_limit):

    '''
    Generates excel sheet for use in importing to tracefinder (via template)
    a, b and y +1 fragments +1 charge
    Chemical formulas


    '''
    list_of_peptides_df = pd.read_csv(list_of_peptides, header=None)
    list_of_peptides = list_of_peptides_df[0].tolist()

    ion_list = []
    mz_list = []
    pdict = {}

    for peptide in list_of_peptides:
        tsg = TheoreticalSpectrumGenerator()
        spec1 = MSSpectrum()
        peptide_to_fragment = AASequence.fromString(peptide)
        p = Param()
        p.setValue("add_a_ions", "true")
        p.setValue("add_metainfo", "true")
        tsg.setParameters(p)
        tsg.getSpectrum(spec1, peptide_to_fragment, 1, 1) # charge range 1:1
        pdict[peptide] = []
        seq = AASequence.fromString(peptide)    
        seq_formula = seq.getFormula()
        pdict[peptide].append(str(seq_formula))
        for ion, peak in zip(spec1.getStringDataArrays()[0], spec1):           
            ion_list.append(ion.decode())
            mz_list.append(peak.getMZ())
            if peak.getMZ() > lower_mass_limit:
                pdict[peptide].append(peak.getMZ())
    list_of_peptides_df = pd.DataFrame.from_dict(pdict, orient='index')
  
    ion_list.clear()
    mz_list.clear()
    pdict.clear()

    with pd.ExcelWriter('SMD_fragments.xlsx') as writer:  
        list_of_peptides_df.to_excel(writer, sheet_name='peptides')

        
