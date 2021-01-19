# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:23:28 2020

@author: Kylie Standage-Beier
"""
import requests #for webpage requests
import json #for loading contents
import time #for timing For-loop API requests
import os 
import numpy as np
#########

"""Main Function, Run button initiates this"""
def design_pegs(data_in):
    parameters = data_in
    data_out = []
    rows = []
    
    """DNA SEQUENCES"""    
    cas9_hp = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
    u6_term = "TTTTTTGTTTT"
    
    """GENERATES DNA REVERSE COMPLEMENT"""
    def dna_rev_comp(DNA):
        dna_a = (DNA.upper().replace(' ',''))  #converts input to uppercase, this is a redundancy
        dna_r = (dna_a[::-1])   #reverse of input DNA
        dna_i = (dna_a.replace("A", "W").replace("T", "X").replace("C", "Y").replace("G", "Z")) #intermediate for making complement
        dna_c = (dna_i.replace("W", "T").replace("X", "A").replace("Y", "G").replace("Z", "C")) #complement
        dna_rc = (dna_c[::-1]) #reverse complement
        return dna_rc, dna_c, dna_r, dna_a #output of function
     
    def GC_content(DNA):
        G_count = 0
        C_count = 0
        A_count = 0
        T_count = 0
        null_count = 0
        GC_count = 0
        
        for i in DNA.upper(): #adds up all G, C, A, T
            if i == 'G':
                G_count += 1
            elif i == 'C':
                C_count += 1
            elif i == 'A':
                A_count += 1
            elif i == 'T':
                T_count += 1
            else:
                null_count += 1
        
        GC_count = G_count + C_count
        GC_percent = ((G_count + C_count)/len(DNA))*100
        return GC_count, GC_percent 
        
    def primer_tm_design(sequence): #enter dna sequence, generates primer with specific Tm, input seq should be 40bp or longer
        tm_set = 63 #Sets goal Tm
        tm = 0 #Starting Tm at begining of For-loop
        #dH is deltaH enthalpy dinucleotide dictionary
        dH = {'AA':-7.9,'TT':-7.9,'AT':-7.2,'TA':-7.2,
              'CA':-8.5,'TG':-8.5,'GT':-8.4,'AC':-8.4,
              'CT':-7.8,'AG':-7.8,'GA':-8.2,'TC':-8.2,
              'CG':-10.6,'GC':-9.8,'GG':-8.0,'CC':-8.0}
        #dS is entropy, -kcal/mol dinucleotide dictionary
        dS = {'AA':-22.6,'TT':-22.6,'AT':-20.4,'TA':-21.3,
              'CA':-22.7,'TG':-22.7,'GT':-22.4,'AC':-22.4,
              'CT':-21.6,'AG':-21.6,'GA':-22.2,'TC':-22.2,
              'CG':-27.2,'GC':-24.4,'GG':-19.9,'CC':-19.9}
                    
        adjust_H = 0.1 
        adjust_S = 2.8
                    
        dH_sum = 0 + adjust_H
        dS_sum = 0 + adjust_S

        for base in range(0,len(sequence)-2):
            dinucleotide = sequence[base:base+2]
            dH_sum += dH[dinucleotide]
            dS_sum += dS[dinucleotide]
            R = 1.987 #universal gas constant 
            M = 500 * 10**-9 #common molarity of primers
            tm_int = ((dH_sum)*1000)/(dS_sum+R*np.log(M))-273.15 #-273.15 converts from kelvin
            tm = (tm_int +16.6*np.log(1))-15 #15 is adjustment to roughly match values
            if tm >= tm_set: #elongates primer until Tm is reached
                break
        
        primer_seq = sequence[0:base+2] #this is final primer sequence        
        return primer_seq
    
    def specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end):
        organism_crispr_dict = {'hg38':'crisprAllTargets','mm10':'crisprAllTargets','sacCer3':'crisprTargets','danRer7':'crisprTargets','ce11':'crisprTargets','dm6':'crisprTargets'}        
        organism_key = guide_api_request_out["genome"] #determines genome of returned reference
        crispr_key = organism_crispr_dict[organism_key] #determines cripsr location in dictionary
        list_guides = guide_api_request_out[crispr_key] #retrieves crispr list
        print("Number of guides: " + str(len(list_guides)))
    
        """ Placeholder lists """     
        guide_index_list = []
        single_guide_list = []
        PAM_list = []
        strand_targeted_list = []
        start_position_list = []
        end_position_list = []
        MIT_spec_score_list = []
        Doen_eff_scor_list = []
        relative_guide_start_list = []
        relative_guide_end_list = []
        
        """ Extract info from CRISPR lists """
        for i in range(len(list_guides)):
            single_guide = list_guides[i]['guideSeq'] #finds guides in list of sub-dictionaries
            PAM = list_guides[i]['pam']
            strand_targeted = list_guides[i]['strand']
            start_position = list_guides[i]['chromStart']
            end_position = list_guides[i]['chromEnd']
            relative_guide_start = 420 - (start_position - input_position + 212) #Corrects coordinates for PAM
            relative_guide_end = end_position - input_position + 208 #Corrects coordinates for PAM
    
            try: #retriev scores for entry
                MIT_spec_score = list_guides[i]['scoreDesc'] #guide_score_list[0].split(":")[1]
                Doench_eff_score = list_guides[i]['fusi'].split('%')[0] #guide_score_list[1].split(":")[1].replace('%','').replace(' ','')
                composite_score = int(MIT_spec_score) + int(Doench_eff_score)
                spec_score = "Specificity: " + str(MIT_spec_score) + " Efficiency: " + Doench_eff_score + " Composite:  " + str(composite_score)
                
            except:
                MIT_spec_score = 0
                Doench_eff_score = 0
                composite_score = 0
                spec_score = "Specificity: " + str(MIT_spec_score) + " Efficiency: " + str(Doench_eff_score)
            
            """ Writes lists for sorting """
            guide_index_list.append(i)
            single_guide_list.append(single_guide)
            PAM_list.append(PAM)
            strand_targeted_list.append(strand_targeted)
            start_position_list.append(start_position)
            end_position_list.append(end_position)
            MIT_spec_score_list.append(MIT_spec_score)
            Doen_eff_scor_list.append(Doench_eff_score)
            relative_guide_start_list.append(relative_guide_start)
            relative_guide_end_list.append(relative_guide_end)
            
            """ Display all guide informationm, for diagnostics"""
            if strand_targeted == '+':
                pos_guides = "Guide: " + str(i) + " " + single_guide + " " + PAM + " " + spec_score + " Strand: " + strand_targeted + " Coordinate: " + str(relative_guide_end)
                #print(pos_guides)
            elif strand_targeted == '-':
                neg_guides = "Guide: " + str(i) + " " + single_guide + " " + PAM + " " + spec_score + " Strand: " + strand_targeted + " Coordinate: " + str(relative_guide_start)
                #print(neg_guides)
    
        """generates list of guides that fall within range"""
        index_list = [] #positve strand guides
        specific_guides_list = []
        specific_guide_end_list = []
        specificity_list = []
        efficiency_list = []
    
        """Generates lists of guides in defined range, generates + and - strand lists"""
        #determines what coordinates to use based off of strand
        if strand == '+':
            coord = relative_guide_end_list
        if strand == '-':
            coord = relative_guide_start_list
        
        for i in range(len(single_guide_list)): #205, 205+((edit_len-2)-(length_of_edit-1))
            if strand_targeted_list[i] == strand: #controls strand of preferred guide
                if int(coord_start) <= coord[i] <= int(coord_end): #numbers will be adapted to match RT, uses direct (non reversed coordinates)
                    index_list.append(i) #makes list of guide indexes that fall within range
                    specific_guides_list.append(single_guide_list[i]) #generates guide sequence list
                    specific_guide_end_list.append(int(coord[i])) #generates end coords
                    specificity_list.append(int(MIT_spec_score_list[i])) #specificity scores for valid guides
                    efficiency_list.append(int(Doen_eff_scor_list[i])) #efficiency scores for valid guides
    
        """Sorts lists of guides within range"""
        """CHANGED!, removed specific_guides_list, and guide_most_specific, int of MIT spec score"""
        ###############
        specificity_list, efficiency_list, specific_guide_end_list, index_list = zip(*sorted(zip(specificity_list, efficiency_list, specific_guide_end_list, index_list))) #returns two or more lists, sorted by first list to last
        ############## 
        guide_index_most_specific = index_list[-1]
        specific_PAM_coord = specific_guide_end_list[-1]

        #return MIT_spec_score_list, relative_guide_start_list, relative_guide_end_list, strand_targeted_list, Doen_eff_scor_list, single_guide_list
        return guide_index_most_specific, specific_PAM_coord

    """DESIGNS POSITIVE STRAND PEGRNA"""
    def pegrna_top(DNA, Edit_site, Edit_nucl, edit_len_a, pbs_len_a, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg, organism, strategy_selected):
        
        """Determines if user is requesting deletion, defines editing product"""
        if "D" in Edit_nucl: #Determines if user is requesting deletion
            deletion_length = int(Edit_nucl[1:len(Edit_nucl)])
            DNA_l = DNA[0:Edit_site]
            DNA_r = DNA[Edit_site+deletion_length:len(DNA)]  
            DNA_edited = DNA_l + DNA_r
        elif "I" in Edit_nucl:
            integraton = str(Edit_nucl[1:len(Edit_nucl)])
            DNA_l = DNA[0:Edit_site+1]
            DNA_r = DNA[1+Edit_site:len(DNA)]
            DNA_edited = DNA_l + integraton + DNA_r
        else: #if deletion not requested, defines edited DNA below
            DNA_l = DNA[0:Edit_site]
            DNA_r = DNA[len(Edit_nucl)+Edit_site:len(DNA)] #Originally DNA[1+Edit_site:len(DNA)], now enables replacement
            DNA_edited = DNA_l + Edit_nucl + DNA_r

        length_of_edit = len(Edit_nucl) #computes length of Edit in function, used to modify PAM search and avoid distant PAMs
        
        DNA_wt_print = DNA[Edit_site-30:Edit_site+30]
        DNA_edit_print = DNA_edited[Edit_site-30:Edit_site+(length_of_edit-1)+30]
        
        """Determines Guide search, specific or near"""
        if design_alg == 'specific':
            strand = '+'
            coord_start = len(DNA)-(205+(edit_len_a-2)-(length_of_edit-1)) #top peg (213-), bot peg (204), edit is 209 
            coord_end = 214 #top peg (213), bot peg (204+), edit it 209
            guide_index_most_specific, specific_PAM_coord = specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end)
            PAM_site = specific_PAM_coord
        else:
            dna_reverse = DNA[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
            pam_rev_search = dna_reverse.index("GG", 205, (205+(edit_len_a-2)-(length_of_edit-1))) 
            #this is reverse search, so end-to-start
            PAM_site = len(DNA)-pam_rev_search-2 #converts back to non-inverted coordinate
                
        protosp_pe = DNA[PAM_site-21:PAM_site-1]
        guide_1 = "G" + DNA[(PAM_site-21):(PAM_site-1)]
        
        ###PBS Length, If blank input, analyzes GC content to determine PBS length###
        if pbs_len_a == "":
            #test 9nt PBS
            pbs9 = DNA[(PAM_site-4-9):(PAM_site-4)]
            GC_count, GC_percent = GC_content(pbs9)
            GC_percent_9 = GC_percent
            #test 12nt PBS
            pbs13 = DNA[(PAM_site-4-13):(PAM_site-4)]
            GC_count, GC_percent = GC_content(pbs13)
            GC_percent_13 = GC_percent
            #test 15nt PBS
            pbs15 = DNA[(PAM_site-4-15):(PAM_site-4)]
            GC_count, GC_percent = GC_content(pbs15)
            GC_percent_15 = GC_percent
            
            if GC_percent_9 >= 60:
                #print("GC > 60")
                pbs_len_a = 9
            elif 40 < GC_percent_13 < 60:
                #print("40 < GC < 60")
                pbs_len_a = 13
            elif GC_percent_15 <= 40:
                #print("GC < 40")
                pbs_len_a = 15
            else:
                pbs_len_a = 13
            pbs_1 = DNA[(PAM_site-4-pbs_len_a):(PAM_site-4)]
        else:
            pbs_1 = DNA[(PAM_site-4-pbs_len_a):(PAM_site-4)]
        
        pbs_top_final = pbs_len_a #to return PBS length
        
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pbs_1)
        pbs_finish = dna_rc
        
        edit_1 = DNA_edited[(PAM_site-4):(PAM_site-4+edit_len_a)] #Generates RT_edit from Edited product
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(edit_1) #Reverse complements RT Edit
        edit_finish = dna_rc #Finished RT_template sequence 
        
        ### pegRNA Oligos ###
        peg_guide_top = "CACC" + guide_1
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(guide_1)
        peg_guide_bot = "AAC" + dna_rc
        peg_edit_top = "TGC" + edit_finish + pbs_finish + u6_term + "CTGCA"
        peg_e_bot_start = edit_finish + pbs_finish + u6_term
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(peg_e_bot_start)
        peg_edit_bot = "G" + dna_rc        
        pegrna_top_finish = str(guide_1 + cas9_hp + edit_finish + pbs_finish + u6_term) #finished pegRNA sequence

        #### Editing strategy function ###
        if strategy_selected == 'PE2':
            pe3_cut_distance = "" #blank since PE2 selected
            pe3_pe_finish = "" #blank since PE2 selected
            pe3_guide = "" #blank since PE2 selected
            ### Oligos ###
            pe3_guide_top = "" #blank since PE2 selected
            pe3_guide_bot = "" #blank since PE2 selected
            print('PE2 running')
            
        elif strategy_selected == 'PE3':
            try: #attempts to find PE3 protospacer
                if design_alg == 'specific':
                    try: #first searches for guides downstream of the edit
                        strand = '-'
                        coord_start = len(DNA)-(PAM_site+90)-2 #
                        coord_end = len(DNA)-(PAM_site+40)-2 #
                        guide_index_most_specific, specific_PAM_coord = specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end)
                        pe3_pam_site = len(DNA)-specific_PAM_coord-2 #converts to correct coordinates
                        print("PE3 specific guide index: " + str(guide_index_most_specific))
                    except: #searches upstream if none found downstream
                        strand = '-'
                        coord_start = len(DNA)-(PAM_site+40)-2 #
                        coord_end = len(DNA)-(PAM_site-90)-2 #
                        guide_index_most_specific, specific_PAM_coord = specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end)
                        pe3_pam_site = len(DNA)-specific_PAM_coord-2 #converts to correct coordinates
                        print("PE3 specific guide index: " + str(guide_index_most_specific))                        
                else:
                    try:
                        pe3_pam_site = DNA.index("CC", PAM_site+40,PAM_site+90) #default +40,+90#searches for PE3 guide PAM on opposite strand
                    except: #searches upstream if none found downstream
                        pe3_pam_site = DNA.index("CC", PAM_site-90,PAM_site-40) 
                        
                pe3_pe_start = DNA[int(pe3_pam_site)+3:int(pe3_pam_site)+23] 
                pe3_cut_distance = int(pe3_pam_site)-int(PAM_site)+10 #
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_pe_start) 
                pe3_pe_finish = dna_rc #
                pe3_guide = "G" + pe3_pe_finish
                ### Oligos ###
                pe3_guide_top = "CACC" + pe3_guide
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_guide)
                pe3_guide_bot = "AAAC" + dna_rc
                print("PE3 running")
            except: #if PE3 protospacer not found
                pe3_cut_distance = ""
                pe3_pe_finish = "PE3 Protospacer not found"
                pe3_guide = ""
                pe3_guide_top = ""
                pe3_guide_bot = ""

        elif strategy_selected == 'PE3B':
            try: #attempts to find PE3B protospacer
                dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA_edited)
                DNA_edited_reverse = dna_r
                pe3_pam_site_init = DNA_edited_reverse.index("CC", Edit_site+3,Edit_site+23) 
                pe3_pam_site = len(DNA_edited)-pe3_pam_site_init-2
                pe3_pe_start = DNA_edited[int(pe3_pam_site)+3:int(pe3_pam_site)+23] 
                pe3_cut_distance = int(pe3_pam_site)-int(PAM_site)+10 #
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_pe_start) 
                pe3_pe_finish = dna_rc #
                pe3_guide = "G" + pe3_pe_finish
                ### Oligos ###
                pe3_guide_top = "CACC" + pe3_guide
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_guide)
                pe3_guide_bot = "AAAC" + dna_rc
            except: #if PE3B protospacer not found
                pe3_cut_distance = "" 
                pe3_pe_finish = "PE3B Protospacer not found"
                pe3_guide = ""
                pe3_guide_top = ""
                pe3_guide_bot = ""
        
        return protosp_pe, PAM_site, guide_1, cas9_hp, edit_finish, pbs_finish, u6_term, pegrna_top_finish, pe3_cut_distance, pe3_pe_finish, pe3_guide, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot, pbs_top_final, DNA_wt_print, DNA_edit_print     
    
    """DESIGNS BOTTOM STRAND pegRNA"""
    def pegrna_bot(DNA, Edit_site, Edit_nucl, edit_len_b, pbs_len_b, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg, organism, strategy_selected):
        
        if "D" in Edit_nucl: #Determines if user is requesting deletion
            deletion_length = int(Edit_nucl[1:len(Edit_nucl)])
            DNA_l = DNA[0:Edit_site]
            DNA_r = DNA[Edit_site+deletion_length:len(DNA)]  
            DNA_edited = DNA_l + DNA_r
        elif "I" in Edit_nucl:
            integraton = str(Edit_nucl[1:len(Edit_nucl)])
            DNA_l = DNA[0:Edit_site+1]
            DNA_r = DNA[1+Edit_site:len(DNA)]
            DNA_edited = DNA_l + integraton + DNA_r
        else: #if deletion not requested, defines edited DNA below
            DNA_l = DNA[0:Edit_site]
            DNA_r = DNA[len(Edit_nucl)+Edit_site:len(DNA)] #Originally DNA[1+Edit_site:len(DNA)], now enables replacement
            DNA_edited = DNA_l + Edit_nucl + DNA_r
        
        length_of_edit = len(Edit_nucl) #computes length of Edit in function, used to modify PAM search and avoid distant PAMs
        
        DNA_wt_print = DNA[Edit_site-30:Edit_site+30]
        DNA_edit_print = DNA_edited[Edit_site-30:Edit_site+(length_of_edit-1)+30]
        
        """Determines PAM coordinates search, specific or near"""
        if design_alg == 'specific':
            dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA)
            dna_rcs = str(dna_rc)
            strand = '-'
            
            #coord_start = len(DNA)-(204+(edit_len_b-2)-(length_of_edit-1)) #top peg (213-), bot peg (204), edit is 209 
            #coord_end = 214 #top peg (213), bot peg (204+), edit it 209     
            """CHANGED!"""
            #######################
            if "D" in Edit_nucl:                
                coord_start = len(DNA)-(204+(edit_len_b-2)-(length_of_edit-1)+deletion_length-1) #top peg (213-), bot peg (204), edit is 209                 
                coord_end = len(DNA)-(204+deletion_length-1)-2    
            elif "I" in Edit_nucl:                        
                coord_start = len(DNA)-(204+(edit_len_b-2)-(length_of_edit-1)) #top peg (213-), bot peg (204), edit is 209 
                coord_end = 214 #top peg (213), bot peg (204+), edit it 209
            else:            
                coord_start = len(DNA)-(204+(edit_len_b-2)-(length_of_edit-1)) #top peg (213-), bot peg (204), edit is 209 
                coord_end = len(DNA)-(204+(length_of_edit-1))-2 #top peg (213), bot peg (204+), edit it 209
            ####################### 
            
            guide_index_most_specific, specific_PAM_coord = specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end)
            print("Bot pegRNA index: " + str(guide_index_most_specific))
            PAM_site_b = specific_PAM_coord
        else:
            dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA)
            dna_rcs = str(dna_rc)
            dna_reverse = dna_rcs[::-1]
            
            """CHANGED!"""
            #######################
            if "D" in Edit_nucl:
                pam_rev_search = dna_reverse.index("GG", (204+deletion_length-1), (204+(edit_len_b-2)-(length_of_edit-1)+deletion_length-1))
            elif "I" in Edit_nucl:
                pam_rev_search = dna_reverse.index("GG", 204, (204+(edit_len_b-2)-(length_of_edit-1)))
            else:
                pam_rev_search = dna_reverse.index("GG", (204+(length_of_edit-1)), (204+(edit_len_b-2)-(length_of_edit-1)))
            #######################            
            PAM_site_b = len(DNA)-pam_rev_search-2 #converts to absolute coordinates, for the reverse complement         
        
        """Designs pegRNA following PAM coordinate"""
        protosp_pe_b = dna_rcs[PAM_site_b-21:PAM_site_b-1]
        guide_1_b = "G" + dna_rcs[(PAM_site_b-21):(PAM_site_b-1)]
        
        ###PBS Length###
        if pbs_len_b == "": #if PBS length is blank, identifies PBS length by GC content
            #test 9nt PBS
            pbs9 = dna_rcs[(PAM_site_b-4-9):(PAM_site_b-4)]
            GC_count, GC_percent = GC_content(pbs9)
            GC_percent_9 = GC_percent
            #test 12nt PBS
            pbs13 = dna_rcs[(PAM_site_b-4-13):(PAM_site_b-4)]
            GC_count, GC_percent = GC_content(pbs13)
            GC_percent_13 = GC_percent
            #test 15nt PBS
            pbs15 = dna_rcs[(PAM_site_b-4-15):(PAM_site_b-4)]
            GC_count, GC_percent = GC_content(pbs15)
            GC_percent_15 = GC_percent
            
            if GC_percent_9 >= 60:
                #print("GC > 60")
                pbs_len_b = 9
            elif 40 < GC_percent_13 < 60:
                #print("40 < GC < 60")
                pbs_len_b = 13
            elif GC_percent_15 <= 40:
                #print("GC < 40")
                pbs_len_b = 15
            else:
                pbs_len_b = 13
            pbs_1 = dna_rcs[(PAM_site_b-4-pbs_len_b):(PAM_site_b-4)] #defining PBS
        else:
            pbs_1 = dna_rcs[(PAM_site_b-4-pbs_len_b):(PAM_site_b-4)] #defining PBS
        
        pbs_bot_final = pbs_len_b
        #pbs_1 = dna_rcs[(PAM_site_b-4-pbs_len_a):(PAM_site_b-4)] #defining PBS
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pbs_1)
        pbs_finish_b = dna_rc
        
        dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA_edited) #defining RT edit
        edit_1 = dna_rc[(PAM_site_b-4):(PAM_site_b-4+edit_len_b)]
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(edit_1)
        edit_finish_b = dna_rc        
        pegrna_bot_finish = str(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term) #full pegRNA seq
        
        ### Oligos ###
        peg_guide_top_b = "CACC" + str(guide_1_b) #Guide oligo top
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(guide_1_b)
        peg_guide_bot_b = "AAC" + dna_rc #Guide oligo bottom
        peg_edit_top_b = "TGC" + edit_finish_b + pbs_finish_b + u6_term + "CTGCA" #Edit oligo top
        peg_e_bot_start = edit_finish_b + pbs_finish_b + u6_term
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(peg_e_bot_start)
        peg_edit_bot_b = "G" + dna_rc #Edit oligo bottom                
        if strategy_selected == 'PE2':
            pe3_cut_distance_b = "" #blank since PE2 selected
            pe3_pe_finish_b = "" #blank since PE2 selected
            pe3_guide_b = "" #blank since PE2 selected
            ### Oligos ###
            pe3_guide_top_b = "" #blank since PE2 selected
            pe3_guide_bot_b = "" #blank since PE2 selected
            
        elif strategy_selected == 'PE3':
            try: #attempts to find PE3 protospacer
                if design_alg == 'specific':
                    try: #first searches for PE3 guide downstream
                        strand = '+'
                        coord_start = len(DNA)-(PAM_site_b)-2-90 #converts pam_site_b into + coordinates range for PE3
                        coord_end = len(DNA)-(PAM_site_b)-2-40 #converts pam_site_b into + coordinates range for PE3
                        guide_index_most_specific, specific_PAM_coord = specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end)
                        pe3_pam_site_b = 420-specific_PAM_coord-2 #converts coordinates match correct
                        print("PE3 Guide index: " + str(guide_index_most_specific))
                    except: #searches upstream if no guide found
                        strand = '+'
                        coord_start = len(DNA)-(PAM_site_b)-2+40 #converts pam_site_b into + coordinates range for PE3
                        coord_end = len(DNA)-(PAM_site_b)-2+90 #converts pam_site_b into + coordinates range for PE3
                        guide_index_most_specific, specific_PAM_coord = specific_guides(guide_api_request_out, input_position, strand, coord_start, coord_end)
                        pe3_pam_site_b = 420-specific_PAM_coord-2 #converts coordinates match correct
                        print("PE3 Guide index: " + str(guide_index_most_specific))
                else:
                    try: #first searches for PE3 guide downstream
                        pe3_pam_site_b = dna_rcs.index("CC", PAM_site_b+40,PAM_site_b+90) #default +40,+90#searches for PE3 guide PAM on opposite strand
                    except: #first searches for PE3 guide downstream
                        pe3_pam_site_b = dna_rcs.index("CC", PAM_site_b-90,PAM_site_b-40) #default +40,+90#searches for PE3 guide PAM on opposite strand
                
                pe3_pe_start_b = dna_rcs[int(pe3_pam_site_b)+3:int(pe3_pam_site_b)+23] #double check spacings
                pe3_cut_distance_b = int(pe3_pam_site_b)-int(PAM_site_b)+10 #
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_pe_start_b) 
                pe3_pe_finish_b = dna_rc #
                pe3_guide_b = "G" + pe3_pe_finish_b
                ### Oligos ###
                pe3_guide_top_b = "CACC" + pe3_guide_b
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_guide)
                pe3_guide_bot_b = "AAAC" + dna_rc
            except: #if PE3 protospacer not found
                pe3_cut_distance_b = ""
                pe3_pe_finish_b = "PE3 Protospacer not found"
                pe3_guide_b = ""
                pe3_guide_top_b = ""
                pe3_guide_bot_b = ""

        elif strategy_selected == 'PE3B':
            try: #attempts to find PE3B protospacer
                dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA_edited)
                DNA_edited_reverse = dna_r
                pe3_pam_site_init = DNA_edited.index("GG", Edit_site,Edit_site+20)
                pe3_pam_site_b = pe3_pam_site_init
                pe3_pe_finish_b = DNA_edited[int(pe3_pam_site_b)-21:int(pe3_pam_site_b)-1] #double check spacings
                pe3_cut_distance_b = int(pe3_pam_site_b)-int(PAM_site_b)+10 #
                pe3_guide_b = "G" + pe3_pe_finish
                ### Oligos ###
                pe3_guide_top_b = "CACC" + pe3_guide_b
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_guide_b)
                pe3_guide_bot_b = "AAAC" + dna_rc
            except: #if PE3B protospacer not found
                pe3_cut_distance_b = "" 
                pe3_pe_finish_b = "PE3B Protospacer not found"
                pe3_guide_b = ""
                pe3_guide_top_b = ""
                pe3_guide_bot_b = ""        
    
        return protosp_pe_b, PAM_site_b, guide_1_b, cas9_hp, edit_finish_b, pbs_finish_b, u6_term, pegrna_bot_finish, pe3_cut_distance_b, pe3_pe_finish_b, pe3_guide_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b, pbs_bot_final, DNA_wt_print, DNA_edit_print   

    """LOOP FOR GENERATING PEGRNA CSV FILE"""
    for index in range(len(parameters)):    
        ### General design is identified on a per-peg basis ##
        organism = parameters[index]['organism'] #'null','string','text','plasmid','hg38','sacCer3','mm10','danRer7','ce11','dm6'
        design_alg = parameters[index]['design_alg'] #'near','specific'
        strategy_selected = parameters[index]['edit_strat'] #'PE2','PE3','PE3B'
        ### Guide parameters ###
        pegRNA_name = str(parameters[index]['name']) #Name of guide in input file    
        chrom_num = str(parameters[index]['chromosome'].replace("Chr","").replace(" ","").replace(",","")) #chromosome of target locus    
        input_position = int(parameters[index]['position'].replace(" ","").replace(",","")) #specific Bp address of locus    
        Edit_nucl = str(parameters[index]['edit'].replace(" ","").replace(",","")).upper() #Edit Specified in RT edit        
        length_of_edit = len(Edit_nucl) #used for modifying PAM site search calculation, avoid spurious ID of distant PAMs        
        ### Statement for RT lengths ###        
        edit_len = parameters[index]['rt'].replace(" ","").replace(",","") #removes spaces and , from input        
        if edit_len == "":
            edit_len_a = edit_len #sets edit_len_a to blank
        else:
            try:
                edit_len_a = int(edit_len) #explicitly attempts to set RT length as an integer
            except:
                edit_len_a = "" #if failiure, sets RT length to blank
        edit_len_b = edit_len_a #sets b bottom to match here. this is to do alternate RT length searches for bottom strand
        ### Statement for PBS lengths ###
        pbs_len = parameters[index]['pbs'].replace(" ","").replace(",","") #Specified Primer length
        if pbs_len == "":
            pbs_len_a = pbs_len #explicitly kepps PBS length as an integer
        else:
            try:
                pbs_len_a = int(pbs_len) #explicitly attempts to set PBS length as an integer
            except:
                pbs_len_a = "" #if failiure, sets PBS length to blank
        pbs_len_b = pbs_len_a
        notes_in = str(parameters[index]['notes']) #This is retained between input and output for note keeping
            
        """Add If statement to change design alg to near if Plasmid or text are used"""
        try:
            ##################### Organism selection #####################
            if organism == "null":
                print("no organism selected")
                DNA = "No Organism, No DNA"
                print(DNA)
            elif organism == "text": #manual 
                file_name = open(parameters[index]['file_dna'], "r")
                print(file_name.read)
                print(file_name.readable()) #TRUE OR FALSE output
                DNA_in = str(file_name.readline()).upper().replace(' ','').replace('\n','') #converts all to upper and removes spaces
                position_start = int(input_position-210) #downstream of edit
                position_end = int(input_position+210) #upstream of edit
                Edit_site = int(209)
                DNA = DNA_in[position_start:position_end]
                design_alg = 'near' #Resets design alg to near if text file used
            elif organism == "plasmid": #Plasmid performs text wrap around if necessary
                file_name = open(parameters[index]['file_dna'], "r")
                print(file_name.read)
                print(file_name.readable()) #TRUE OR FALSE output
                DNA_in = str(file_name.readline()).upper().replace(' ','').replace('\n','') #converts all to upper and removes spaces and new lines
                Edit_site = int(209)
                design_alg = 'near' #Resets design alg to near if text file used
                #determines position in plasmid, rotates plasmid to accomodate edits
                if 210 <= input_position <= len(DNA_in)-210:
                    position_start = int(input_position-210) #downstream of edit
                    position_end = int(input_position+210) #upstream of edit
                    DNA = DNA_in[position_start:position_end]
                elif input_position < 210:
                    plasmid_shift = str(DNA_in[len(DNA_in)-500:len(DNA_in)] + DNA_in[0:len(DNA_in)-500])
                    position_start = int(input_position-210+500) #downstream of edit
                    position_end = int(input_position+210+500) #upstream of edit
                    DNA = plasmid_shift[position_start:position_end]
                elif len(DNA_in)-input_position < 210:
                    plasmid_shift = str(DNA_in[500:len(DNA_in)] + DNA_in[0:500])
                    position_start = int(input_position-210-500) #downstream of edit
                    position_end = int(input_position+210-500) #upstream of edit
                    DNA = plasmid_shift[position_start:position_end]
                elif organism == 'string':
                    DNA_in = str(parameters[index]['file_dna']).upper().replace(' ','').replace('\n','') #converts all to upper and removes spaces
                    position_start = int(input_position-210) #downstream of edit
                    position_end = int(input_position+210) #upstream of edit
                    Edit_site = int(209)
                    DNA = DNA_in[position_start:position_end]
                    design_alg = 'near' #Resets design alg to near if text file used
                
            else: #Human Homo sapiens Hg38 and other
                roman_numeral = {'1':'I','2':'II','3':'III','4':'IV','5':'V','6':'VI','7':'VII','8':'VIII','9':'IX','10':'X','11':'XI','12':'XII','13':'XIII','14':'XIV','15':'XV','16':'XVI','M':'M'} #M is for mitochondira
                if organism == 'sacCer3': #converts numerical input into roman numerals for API request
                    chrom_num_req = roman_numeral[chrom_num]
                elif organism == 'ce11':
                    chrom_num_req = roman_numeral[chrom_num]
                else:
                    chrom_num_req = chrom_num
                        
                ### API Request for reference genome ###
                position_start = str(input_position-210) #downstream of edit
                position_end = str(input_position+210) #upstream of edit
                url_req = str("https://api.genome.ucsc.edu/getData/sequence?genome=" + organism + ";chrom=chr" + chrom_num_req + ";start=" + position_start +";end=" + position_end) #URL of UCSC Genome browaer Hg38 API
                print(url_req)
                api_request = requests.get(url_req) #Begins API request
                api = json.loads(api_request.content)
                DNA_in = api['dna']
                DNA = DNA_in.upper().replace(' ','') #converts all to upper and removes spaces
                Edit_site = int(209) #this is hard coded position of edit in seq returned from API
                time.sleep(0.5)
                print(DNA)
                
            ##################### Test if pegRNA can be designed, Positive Strand #####################
            """Initial PAM search: determines if any valid guides available, Determines RT length if necessary"""
            if edit_len_a == "": #if blank RT length input, searches for correct RT length
                for i in range(23):
                    try:
                        RT_test = i+7 #will test RT lengths from 7 to 23
                        dna_reverse = DNA[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                        pam_rev_search = dna_reverse.index("GG", 205, 205+((int(RT_test)-2)-(length_of_edit-1))) #determines PAM site from reverse, to place as close to edit as possible 
                        PAM_site = len(DNA)-pam_rev_search-2 #converts PAM rev search to absolute coordinates
                        break
                    except:
                        null = "null"
                edit_len_a = RT_test+3 #+3 adds extra pad to RT
            else:
                try:
                    dna_reverse = DNA[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                    pam_rev_search = dna_reverse.index("GG", 205, 205+((int(edit_len_a)-2)-(length_of_edit-1))) #determines PAM site from reverse, to place as close to edit as possible 
                    PAM_site = len(DNA)-pam_rev_search-2 #converts PAM rev search to absolute coordinates
                except:
                    print("PAM not found (Check RT setting)")
                
            """Guide specificity search: retrieves specific guides flanking edit if selected"""
            if design_alg == 'specific':
                guide_retrieve_start = input_position - 119 #downstream bound of guide retrieval
                guide_retrieve_end = input_position + 120 #upstream bound of guide retrieval
                organism_crispr_dict = {'hg38':'crisprAllTargets','mm10':'crisprAllTargets','sacCer3':'crisprTargets','danRer7':'crisprTargets','ce11':'crisprTargets','dm6':'crisprTargets'}         
                guide_url_req = str("https://api.genome.ucsc.edu/getData/track?genome=" + organism + ";track=" + str(organism_crispr_dict[organism]) + ";chrom=chr" + chrom_num + ";start=" + str(guide_retrieve_start) + ";end=" + str(guide_retrieve_end)) #URL of UCSC Genome browaer Hg38 API
                #print(guide_url_req)
                guide_api_request = requests.get(guide_url_req)
                guide_api_request_out = json.loads(guide_api_request.content)
                time.sleep(0.5)
                #print(guide_api_request_out[organism_crispr_dict[organism]])
            else:
                guide_api_request_out = ''
                
            """Design top pegRNA function"""
            protosp_pe, PAM_site, guide_1, cas9_hp, edit_finish, pbs_finish, u6_term, pegrna_top_finish, pe3_cut_distance, pe3_pe_finish, pe3_guide, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot, pbs_top_final, DNA_wt_print, DNA_edit_print = pegrna_top(DNA, Edit_site, Edit_nucl, edit_len_a, pbs_len_a, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg, organism, strategy_selected)
            
            print("Top strand protospacer " + str(index) + ":" + protosp_pe)
            print("Top pEGRNA for " + str(index) + ":")
            print("pegRNA for "  + str(index) + ":")
            #pegrna_top = str(guide_1 + cas9_hp + edit_finish + pbs_finish + u6_term)        
            print(guide_1 + cas9_hp + edit_finish + pbs_finish + u6_term) #print pegRNA for troubleshooting
            pegrna_top_finish = str(guide_1 + cas9_hp + edit_finish + pbs_finish + u6_term) #explicitly places parts together to make pegRNA
                
                
            ### Primers ###
            sequence = str(DNA[0:60])
            primer_seq = primer_tm_design(sequence) #forward primer, tm set
            seq_f_finish = primer_seq
            
            dna_rc , dna_c, dna_r, dna_a  = dna_rev_comp(DNA[len(DNA)-60:len(DNA)]) #reverse primer, tm set
            sequence = str(dna_rc)
            primer_seq = primer_tm_design(sequence)
            seq_r_finish = primer_seq
            
            ### Edit Defined ###
            if "D" in Edit_nucl: #Determines if user is requesting deletion
                editing_product = str(Edit_nucl)
            elif "I" in Edit_nucl:
                editing_product = str(Edit_nucl)
            else:
                editing_product = str(DNA[209:210+(len(Edit_nucl)-1)]) + "-to-" + str(Edit_nucl) #generates X-to-Y format for edit annotation
                
            print("PE3 Cut Distance: " + str(pe3_cut_distance))
            print("PE3 Protospacer: " + pe3_pe_finish)
            print("Forward Sequencing: " + seq_f_finish)
            print("Reverse Sequencing: " + seq_r_finish)
            print("Length of Edit " + str(index) + ": " + str(length_of_edit))
            strand_t = "Positive"
            #validation_top = "valid"
            try: #tries to find poly-T track
                poly_t_check =  str(edit_finish + pbs_finish).index("TTTT")
                notes = notes_in + ' Contains Poly-U(T) track' #displays warning for poly-T track
            except:
                notes = notes_in
            
        except ValueError:
            pegrna_top_finish = "No Valid pegRNA"
            strand_t = "Positive"
            seq_f_finish = " "
            seq_r_finish = " "
            editing_product = " "
            protosp_pe = " "
            edit_finish = " "
            pbs_top_final = " "
            pbs_finish = " "
            pe3_pe_finish = " "
            pe3_cut_distance = " "
            peg_guide_top = " "
            peg_guide_bot = " "
            peg_edit_top = " "
            peg_edit_bot = " "
            pe3_guide_top = " "
            pe3_guide_bot = " "
            notes = notes_in
            DNA_wt_print = " "
            DNA_edit_print = " "
            print("No Valid positive strand pEGRNAs for " + str(index))
        
        row_t = {'index':index,'name':pegRNA_name, 'chromosome':chrom_num, 'position':input_position, 'product':editing_product, 'strand':strand_t, 'pegrna':pegrna_top_finish, 'proto_peg':protosp_pe, 'rt':edit_finish, 'rt_length':edit_len_a, 'pbs':pbs_finish, 'pbs_length':pbs_top_final, 'proto_pe3':pe3_pe_finish, 'pe3_distance':pe3_cut_distance, 'peg_guide_top':peg_guide_top, 'peg_guide_bot':peg_guide_bot, 'peg_ext_top':peg_edit_top, 'peg_ext_bot':peg_edit_bot, 'pe3_top':pe3_guide_top, 'pe3_bot':pe3_guide_bot, 'primer_f':seq_f_finish, 'primer_r':seq_r_finish, 'notes':notes, 'wt_seq':DNA_wt_print, 'edit_seq':DNA_edit_print}
        
        rows.append(row_t)
        #Negative strand pegRNAs below
        try:
            """Identifies if PAM in viable position"""
            dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA)
            dna_rcs = str(dna_rc)
            """CHANGED!"""
            #######################
            ### Identifies RT length for bottom strand if undefined ###
            if edit_len_b == "": #if blank RT length input, searches for correct RT length
                for i in range(23):
                    try:
                        RT_test = i+7 #will test RT lengths from 10 to 30
                        dna_reverse = dna_rcs[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                        if "D" in Edit_nucl:
                            deletion_length = int(Edit_nucl[1:len(Edit_nucl)])
                            pam_rev_search = dna_reverse.index("GG", (204+deletion_length-1), (204+(RT_test-2)-(length_of_edit-1))) #determines if PAM is within 
                        elif "I" in Edit_nucl:
                            pam_rev_search = dna_reverse.index("GG", 204, (204+(RT_test-2)-(length_of_edit-1)))
                        else:
                            pam_rev_search = dna_reverse.index("GG", (204+(length_of_edit-1)), (204+(RT_test-2)-(length_of_edit-1)))
                        
                        PAM_site_b = len(DNA)-pam_rev_search-2
                        break
                    except:
                        null = "null"
                
                edit_len_b = RT_test+3 #adds extra pad to RT
            else:
                try:
                    dna_reverse = dna_rcs[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                    if "D" in Edit_nucl:
                        deletion_length = int(Edit_nucl[1:len(Edit_nucl)])
                        pam_rev_search = dna_reverse.index("GG", (204+deletion_length), (204+(edit_len_b-2)-(length_of_edit-1)+deletion_length)) #determines if PAM is within 
                    elif "I" in Edit_nucl:
                        pam_rev_search = dna_reverse.index("GG", 204, (204+(edit_len_b-2)-(length_of_edit-1)))
                    else:
                        pam_rev_search = dna_reverse.index("GG", (204+(length_of_edit-1)), (204+(edit_len_b-2)-(length_of_edit-1)))
                        #pam_rev_search = dna_reverse.index("GG", 204, (204+(int(edit_len_b)-2)-(length_of_edit-1))) #determines if PAM is within                         
                    PAM_site_b = len(DNA)-pam_rev_search-2
                except:
                    print("PAM not found (Check RT setting)")
                #######################
            """Design bottom pegRNA function"""
            protosp_pe_b, PAM_site_b, guide_1_b, cas9_hp, edit_finish_b, pbs_finish_b, u6_term, pegrna_bot_finish, pe3_cut_distance_b, pe3_pe_finish_b, pe3_guide_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b, pbs_bot_final, DNA_wt_print, DNA_edit_print = pegrna_bot(DNA, Edit_site, Edit_nucl, edit_len_b, pbs_len_b, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg, organism, strategy_selected)
            
            print("Bottom strand protospacer " + str(index) + ":" + protosp_pe_b)
            print("Bottom pEGRNA for " + str(index) + ":")
            print("pegRNA for "  + str(index) + ":")
            #pegrna_bot_finish = str(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term)
            print(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term)
            pegrna_bot_finish = str(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term)
                
            ### Primers ###
            sequence = str(DNA[0:60])
            primer_seq = primer_tm_design(sequence) #forward primer, tm set
            seq_f_finish = primer_seq
                
            dna_rc , dna_c, dna_r, dna_a  = dna_rev_comp(DNA[len(DNA)-60:len(DNA)]) #reverse primer, tm set
            sequence = str(dna_rc)
            print(sequence)
            primer_seq = primer_tm_design(sequence)
            seq_r_finish = primer_seq

            if "D" in Edit_nucl: #Determines if user is requesting deletion
                editing_product = str(Edit_nucl)
            elif "I" in Edit_nucl:
                editing_product = str(Edit_nucl)
            else:
                editing_product = str(DNA[209:210+(len(Edit_nucl)-1)]) + "-to-" + str(Edit_nucl) #generates X-to-Y format for edit annotation
                
            print("PE3 Cut Distance: " + str(pe3_cut_distance_b))
            print("PE3 Protospacer: " + pe3_pe_finish_b)
            print("Forward Sequencing: " + seq_f_finish)
            print("Reverse Sequencing: " + seq_r_finish)
            #validation_bot = "valid"
            strand_b = "Negative"
            try: #tries to find poly-T track
                poly_t_check =  str(edit_finish_b + pbs_finish_b).index("TTTT")
                notes = notes_in + ' Contains Poly-U(T) track' #displays warning for poly-T track
            except:
                notes = notes_in
        except ValueError:
            pegrna_bot_finish = "No Valid pegRNA"
            editing_product = " "
            strand_b = "Negative"
            seq_f_finish = " "
            seq_r_finish = " "
            seq_r_finish = " "
            protosp_pe_b = " "
            edit_finish_b = " "
            pbs_finish_b = " "
            pbs_bot_final = " "
            pe3_pe_finish_b = " "
            pe3_cut_distance_b = " "
            peg_guide_top_b = " "
            peg_guide_bot_b = " "
            peg_edit_top_b = " "
            peg_edit_bot_b = " "
            pe3_guide_top_b = " "
            pe3_guide_bot_b = " "
            notes = notes_in
            DNA_wt_print = " "
            DNA_edit_print = " "
            #validation_bot = "error"
            print("No Valid negative strand  pEGRNAs " + str(index))        

        row_b = {'index':index, 'name':pegRNA_name, 'chromosome':chrom_num, 'position':input_position, 'product':editing_product, 'strand':strand_b, 'pegrna':pegrna_bot_finish, 'proto_peg':protosp_pe_b, 'rt':edit_finish_b, 'rt_length':edit_len_b, 'pbs':pbs_finish_b, 'pbs_length':pbs_bot_final, 'proto_pe3':pe3_pe_finish_b, 'pe3_distance':pe3_cut_distance_b, 'peg_guide_top':peg_guide_top_b, 'peg_guide_bot':peg_guide_bot_b, 'peg_ext_top':peg_edit_top_b, 'peg_ext_bot':peg_edit_bot_b, 'pe3_top':pe3_guide_top_b, 'pe3_bot':pe3_guide_bot_b, 'primer_f':seq_f_finish, 'primer_r':seq_r_finish, 'notes':notes, 'wt_seq':DNA_wt_print, 'edit_seq':DNA_edit_print}
        rows.append(row_b)
        
    data_out = rows
    
    return data_out
