import itertools
from pyopenms import *
import pandas as pd


def dpl_SMD_using_dipeptides(starting_dipeptides, lower_mass_limit):

    '''
    Generates excel sheet for use in importing to tracefinder (via template)
    a, b and y +1 fragments +1 charge
    Chemical formulas


    '''

    #---Generate tetrapeptides
    tetrapeptide_tlist = list(itertools.product(starting_dipeptides, starting_dipeptides))
    tetrapeptide_list = [''.join(tups) for tups in tetrapeptide_tlist]
    
    #---Generate hexapeptides
    hexapeptide_tlist = list(itertools.product(tetrapeptide_list, starting_dipeptides))
    hexapeptide_list = [''.join(tups) for tups in hexapeptide_tlist]

    #---Generate octapeptides
    octapeptide_tlist = list(itertools.product(tetrapeptide_list, tetrapeptide_list))
    octapeptide_list = [''.join(tups) for tups in octapeptide_tlist]

    #---Generate ZXZ tripeptides
    tripeptide_list = []
    for peptide in tetrapeptide_list:
        tripeptide_list.append(peptide[0:3])

    #---Generate XZX tripeptides
    for peptide in tetrapeptide_list:
        tripeptide_list.append(peptide[1:4])
        
    #---Generate ZXZXZ pentapeptides
    pentapeptide_list = []
    for peptide in hexapeptide_list:
        pentapeptide_list.append(peptide[0:5])

    #---Generate XZXZX pentapeptides
    for peptide in hexapeptide_list:
        pentapeptide_list.append(peptide[1:6])

    #---Generate ZXZXZXZ heptapeptides
    heptapeptide_list = []
    for peptide in octapeptide_list:
        heptapeptide_list.append(peptide[0:6])
        
    #---Generate XZXZXZX heptapeptides
    for peptide in octapeptide_list:
        heptapeptide_list.append(peptide[1:7])

    ion_list = []
    mz_list = []
    pdict = {}

    for peptide in starting_dipeptides:
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
    starting_dipeptides_df = pd.DataFrame.from_dict(pdict, orient='index')
  
    ion_list.clear()
    mz_list.clear()
    pdict.clear()

    for peptide in tripeptide_list:
        tsg = TheoreticalSpectrumGenerator()
        spec1 = MSSpectrum()
        peptide_to_fragment = AASequence.fromString(peptide)
        # standard behavior is adding b- and y-ions of charge 1
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
    tripeptide_df = pd.DataFrame.from_dict(pdict, orient='index')

    #~~~Generate a,b and y fragments and chemical formulas from lists of peptides
    
    ion_list.clear()
    mz_list.clear()
    pdict.clear()

    for peptide in tetrapeptide_list:
        tsg = TheoreticalSpectrumGenerator()
        spec1 = MSSpectrum()
        peptide_to_fragment = AASequence.fromString(peptide)
        # standard behavior is adding b- and y-ions of charge 1
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
    tetrapeptide_df = pd.DataFrame.from_dict(pdict, orient='index')
    
    ion_list.clear()
    mz_list.clear()
    pdict.clear()

    for peptide in pentapeptide_list:
        tsg = TheoreticalSpectrumGenerator()
        spec1 = MSSpectrum()
        peptide_to_fragment = AASequence.fromString(peptide)
        # standard behavior is adding b- and y-ions of charge 1
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
    pentapeptide_df = pd.DataFrame.from_dict(pdict, orient='index')
    
    ion_list.clear()
    mz_list.clear()
    pdict.clear()

    for peptide in hexapeptide_list:
        tsg = TheoreticalSpectrumGenerator()
        spec1 = MSSpectrum()
        peptide_to_fragment = AASequence.fromString(peptide)
        # standard behavior is adding b- and y-ions of charge 1
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
    hexapeptide_df = pd.DataFrame.from_dict(pdict, orient='index')
    
    ion_list.clear()
    mz_list.clear()
    pdict.clear()
    for peptide in heptapeptide_list:
        tsg = TheoreticalSpectrumGenerator()
        spec1 = MSSpectrum()
        peptide_to_fragment = AASequence.fromString(peptide)
        # standard behavior is adding b- and y-ions of charge 1
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
    heptapeptide_df = pd.DataFrame.from_dict(pdict, orient='index')
    
    ion_list.clear()
    mz_list.clear()
    pdict.clear()   
    for peptide in octapeptide_list:
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
    octapeptide_df = pd.DataFrame.from_dict(pdict, orient='index')

    with pd.ExcelWriter('dipeptide_DPL_fragments.xlsx') as writer:  
        starting_dipeptides_df.to_excel(writer, sheet_name='dipeptides')
        tripeptide_df.to_excel(writer, sheet_name='tripeptides')
        tetrapeptide_df.to_excel(writer, sheet_name='tetrapeptides')
        pentapeptide_df.to_excel(writer, sheet_name='pentapeptides')
        hexapeptide_df.to_excel(writer, sheet_name='hexapeptides')
        heptapeptide_df.to_excel(writer, sheet_name='heptapeptides')
        octapeptide_df.to_excel(writer, sheet_name='octapeptides')
        
