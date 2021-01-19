# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:23:28 2020

@author: Kylie Standage-Beier
"""
import csv #for handling input and output CSVs
import requests #for webpage requests
import json #for loading contents
import time #for timing For-loop API requests
from tkinter import * #For UI
from PIL import ImageTk,Image #For displaying images
from tkinter import filedialog #for file dialog interface
from tkinter import ttk #needed for combo boxes
import os 
import sys
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import colors as mcolors
import numpy as np
#########



print("Hello")

root = Tk()
root.title('PINE CONE')
exe_path = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1]) #is able to return correct path to executable

logos_path = os.path.join(exe_path,'logos')
program_logo_path = os.path.join(logos_path,'pineconelogo.icns')
root.iconbitmap(program_logo_path) #logo
root.geometry("520x550") #Size of UI window

"""OPENS CSV INPUT PARAMETERS"""
input_file_path = "aaaaaa"
input_file_path_text = "zzzzz"
parameters = "bbbbbb"
organism = "null"
analysis_button = "null"
strategy_selected = 'PE2_a' #sets as default value for strat button 
design_alg = 'near' #sets as default value for guide design button 

"""UI Stuff"""
"""Organism Dropdown menu, selection dictionary defines organism searched"""
def comboclick(event):
    global label2
    global organism    
    print(myCombo.get())
    
    label2.grid_forget()
    #dictionary converts text into general term, t be used elsewhere (for more complicated programs)
    selction_dictionary = {'Select Organism':'null',
                           'Manual (.txt)':'text',
                           "Plasmid (.txt)":'plasmid',
                           'Human (hg38)':'hg38', 
                           'Yeast (S288C)':'sacCer3',
                           'Mouse (mm10)':'mm10',
                           'Zebrafish (danRer7)':'danRer7',
                           'Roundworm (ce11)':'ce11',
                           'Fruitfly (dm6)':'dm6'}
    
    selection_intermediate = selction_dictionary[myCombo.get()] #just an intermediate to translate between dictionaries
    organism = selction_dictionary[myCombo.get()] #organism is name of reference genome as appears in URL
    #dictionary converts selection variable to image name
    organism_icons = {'null':'Query_logo.png', 
                      'hg38':'Human_logo.png', 
                      'sacCer3':'yeast_logo.png',
                      'mm10':'mouse_logo.png',
                      'text':'manual_text_logo.png',
                      'plasmid':'manual_text_logo.png',
                      'danRer7':'zebrafish_logo.png',
                      'ce11':'roundword_logo.png',
                      'dm6':'fruitfly_logo.png'}
    
    image_to_display = os.path.join(logos_path, organism_icons[selection_intermediate])

    if myCombo.get() == "Manual (.txt)":
        my_btn = Button(root, text="Select DNA Sequence (.txt)", command=select_text_file).grid(row=9, column=0) #if manual selected creates file dialog for text input
    elif myCombo.get() == "Plasmid (.txt)":
        my_btn = Button(root, text="Select DNA Sequence (.txt)", command=select_text_file).grid(row=9, column=0)
    else:
        status_label_organism = "Organism: " + myCombo.get()
    image1 = ImageTk.PhotoImage(Image.open(image_to_display)) #image for organism selection
    label2 = Label(image=image1)
    label2.image = image1
    label2.grid(row=4, column=0, rowspan=2) 

"""Guide Preference Dropdown menu, defines design_alg"""
def comboclick1(event):
    global design_alg
    #Guide prefence dictionary, near is default. also set above
    guide_pref_dict = {'Select Protospacer Preference':'near',
                           'Nearest Protospacers to Edit':'near',
                           "High Specificity (Requires Ref.Genome)":'specific'}
    design_alg = guide_pref_dict[myCombo1.get()]
    print("Guide Design Algorithm: " + design_alg)

"""PE Editing Strategy Radio Buttons Function"""
MODES = [("PE2","PE2_a"),
         ("PE3","PE3_a"),
         ("PE3B","PE3B_a")]

strategy_input = StringVar()
strategy_input.set("PE2_a") #sets default mode


strategy_logo_path = os.path.join(logos_path,'PE2_Logo.png')
image_strategy = ImageTk.PhotoImage(Image.open(strategy_logo_path))
image_label_strategy = Label(image=image_strategy)
image_label_strategy.image = image_strategy
image_label_strategy.grid(row=5, column=1, columnspan=3)

"""Select .CSV file"""
def select_file():
    global input_file_path
    root.filename = filedialog.askopenfilename(initialdir = "/", title="Select File", filetypes = (("csv files", "*.csv"),("all files","*.*")))
    csv_file_name = os.path.basename(root.filename)
    csv_file_display = "Edit File: " + csv_file_name
    my_label = Label(root, text=csv_file_display).grid(row=8, column=0, columnspan=1) #sticky=W
    print (root.filename)
    input_file_path = root.filename

"""If manual input selected, file select for .txt file"""
def select_text_file():
    global input_file_path_text
    root.filename_text = filedialog.askopenfilename(initialdir = "/", title="Select File", filetypes = (("text files", "*.txt"),("all files","*.*")))
    text_file_name = os.path.basename(root.filename_text)
    text_file_display = "DNA Sequence: " + text_file_name
    my_label = Label(root, text=text_file_display).grid(row=10, column=0, columnspan=1)#sticky=W
    input_file_path_text = root.filename_text
    #print(root.filename_text)

"""Main Function, Run button initiates this"""
def run_button():
    global parameters
    global organism
    global output_file_name
    global output_file_name_no_ext
    global design_alg
    
    if output_file_entry.get() == "": #Sets default name of output if no output name is entered
        output_file_name_no_ext = "pinecone_output"
        output_file_name = output_file_name_no_ext + ".csv"
    else: #Sets output name to user input
        output_file_name_no_ext = output_file_entry.get()
        output_file_name = output_file_name_no_ext + ".csv"
        
    output_file_status = "'" + output_file_name + "' is complete"
    clicklabel = Label(root, text=output_file_status, fg='green', font=('helvetica', 12, 'bold'))
    clicklabel.grid(row=8, column=1, columnspan=3)
    retrieve_file_path =  input_file_path
    
    with open(retrieve_file_path) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        data = [list(row) for row in readCSV]
    print(data)
    parameters = data
        
    """DNA SEQUENCE INPUT (old)"""    
    #DNA = input("Enter DNA Sequence: ").replace(" ", "").upper()
    #Edit_site = int(input("Enter Nucleotide Position of Edit (Bp): "))-1
    #Edit_nucl = input("Enter Desired Nucleotide Edit (A,T,C,G): ").upper()
    #edit_len = int(input("Primer Editing Area Length (def. 10 Bp): "))
    #pbs_len = int(input("Primer Binding Site Length (def. 13 Bp): "))
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
        #print("GC Count: " + str(GC_count))
        #print("GC Percent: " + str(GC_percent)+"%")
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
        #print("Primer legnth: "+ str(len(primer_seq)))
        #print("Primer Seq: " + primer_seq)
        #print("Calculated TM is: " + str(tm))
        
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
        #positive strand sort
        """CHANGED!, removed specific_guides_list, and guide_most_specific, int of MIT spec score"""
        ###############
        specificity_list, efficiency_list, specific_guide_end_list, index_list = zip(*sorted(zip(specificity_list, efficiency_list, specific_guide_end_list, index_list))) #returns two or more lists, sorted by first list to last
        ##############          
        guide_index_most_specific = index_list[-1]
        specific_PAM_coord = specific_guide_end_list[-1]
        return guide_index_most_specific, specific_PAM_coord

    """DESIGNS POSITIVE STRAND PEGRNA"""
    def pegrna_top(DNA, Edit_site, Edit_nucl, edit_len_a, pbs_len_a, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg):
        
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
        if strategy_selected == 'PE2_a':
            pe3_cut_distance = "" #blank since PE2 selected
            pe3_pe_finish = "" #blank since PE2 selected
            pe3_guide = "" #blank since PE2 selected
            ### Oligos ###
            pe3_guide_top = "" #blank since PE2 selected
            pe3_guide_bot = "" #blank since PE2 selected
            print('PE2 running')
            
        elif strategy_selected == 'PE3_a':
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

        elif strategy_selected == 'PE3B_a':
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
    def pegrna_bot(DNA, Edit_site, Edit_nucl, edit_len_b, pbs_len_b, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg):
        
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
        if strategy_selected == 'PE2_a':
            pe3_cut_distance_b = "" #blank since PE2 selected
            pe3_pe_finish_b = "" #blank since PE2 selected
            pe3_guide_b = "" #blank since PE2 selected
            ### Oligos ###
            pe3_guide_top_b = "" #blank since PE2 selected
            pe3_guide_bot_b = "" #blank since PE2 selected
            
        elif strategy_selected == 'PE3_a':
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

        elif strategy_selected == 'PE3B_a':
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
    output_file_name_path = exe_path.replace('','') + '/' +output_file_name #app_path sets write file path to same directory as app
    
    with open(output_file_name_path,'w') as f1: #Creates CSV file
        writer=csv.writer(f1, delimiter=',',lineterminator='\n',) #assigns writer parameters, \t
        header = ['Name','Chromosome','Position', 'Edit', 'Strand', 'pegRNA','Protospacer' ,'RT Template','RT Length','Primer Sequence','PBS Length', 'PE3 or 3B Protospacer', 'PE3 or 3B Distance', 'peg_Guide_top', 'peg_Guide_bot', 'peg_Edit_top', 'peg_Edit_bot', 'PE3_Guide_top', 'PE3_Guide_bot', 'Seq_F', 'Seq_R', 'Notes', 'DNA WT (+/-30 Bp)', 'DNA Edited (+/-30 Bp)'] #header of CSV file output
        writer.writerow(header) #creates header in file
        for index in range(1,len(parameters)):
            ###Input_CSV###
            pegRNA_name = str(parameters[index][0]) #Name of guide in input file
            chrom_num = str(parameters[index][1].replace("Chr","").replace("Chr.","").replace(" ","").replace(",","")) #chromosome of target locus
            input_position = int(parameters[index][2].replace(" ","").replace(",","")) #specific Bp address of locus
            Edit_nucl = str(parameters[index][3].replace(" ","").replace(",","")).upper() #Edit Specified in RT edit
            length_of_edit = len(parameters[index][3]) #used for modifying PAM site search calculation, avoid spurious ID of distant PAMs
            
            ### Statement for RT lengths ###
            edit_len = parameters[index][4].replace(" ","").replace(",","") #removes spaces and , from input
            if edit_len == "":
                edit_len_a = edit_len #sets edit_len_a to blank
            else:
                try:
                    edit_len_a = int(edit_len) #explicitly attempts to set RT length as an integer
                except:
                    edit_len_a = "" #if failiure, sets RT length to blank
            
            edit_len_b = edit_len_a #sets b bottom to match here. this is to do alternate RT length searches for bottom strand
            
            ### Statement for PBS lengths ###
            pbs_len = parameters[index][5].replace(" ","").replace(",","") #Specified Primer length
            if pbs_len == "":
                pbs_len_a = pbs_len #explicitly kepps PBS length as an integer
            else:
                try:
                    pbs_len_a = int(pbs_len) #explicitly attempts to set PBS length as an integer
                except:
                    pbs_len_a = "" #if failiure, sets PBS length to blank
            
            pbs_len_b = pbs_len_a
            
            notes_in = str(parameters[index][6]) #This is retained between input and output for note keeping
            
            """Add If statement to change design alg to near if Plasmid or text are used"""
            try:
                ##################### Organism selection #####################
                if organism == "null":
                    print("no organism selected")
                    DNA = "No Organism, No DNA"
                    print(DNA)
                elif organism == "text": #manual 
                    file_name = open(input_file_path_text, "r")
                    print(file_name.read)
                    print(file_name.readable()) #TRUE OR FALSE output
                    DNA_in = str(file_name.readline()).upper().replace(' ','').replace('\n','') #converts all to upper and removes spaces
                    position_start = int(input_position-210) #downstream of edit
                    position_end = int(input_position+210) #upstream of edit
                    Edit_site = int(209)
                    DNA = DNA_in[position_start:position_end]
                    design_alg = 'near' #Resets design alg to near if text file used
                elif organism == "plasmid": #Plasmid performs text wrap around if necessary
                    file_name = open(input_file_path_text, "r")
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
                
                #dna_reverse = DNA[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                #pam_rev_search = dna_reverse.index("GG", 205, 205+((edit_len-2)-(length_of_edit-1))) #determines PAM site from reverse, to place as close to edit as possible 
                #PAM_site = len(DNA)-pam_rev_search-2
                #print(PAM_site)
                
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
                protosp_pe, PAM_site, guide_1, cas9_hp, edit_finish, pbs_finish, u6_term, pegrna_top_finish, pe3_cut_distance, pe3_pe_finish, pe3_guide, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot, pbs_top_final, DNA_wt_print, DNA_edit_print = pegrna_top(DNA, Edit_site, Edit_nucl, edit_len_a, pbs_len_a, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg)
                
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
            row_t = [pegRNA_name, chrom_num, input_position, editing_product, strand_t, pegrna_top_finish, protosp_pe, edit_finish, edit_len_a, pbs_finish, pbs_top_final, pe3_pe_finish, pe3_cut_distance, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot, seq_f_finish, seq_r_finish, notes, DNA_wt_print, DNA_edit_print] 
            writer.writerow(row_t) #writes new row eachtime
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
                                print("Deletion Length is: " + str(deletion_length))
                                pam_rev_search = dna_reverse.index("GG", (204+deletion_length-1), (204+(RT_test-2)-(length_of_edit-1))) #determines if PAM is within 
                                print(dna_reverse[(204+deletion_length-1):(204+(RT_test-2)-(length_of_edit-1)+deletion_length-1)])
                            elif "I" in Edit_nucl:
                                pam_rev_search = dna_reverse.index("GG", 204, (204+(RT_test-2)-(length_of_edit-1)))
                            else:
                                pam_rev_search = dna_reverse.index("GG", (204+(length_of_edit-1)), (204+(RT_test-2)-(length_of_edit-1)))
                            
                            #pam_rev_search = dna_reverse.index("GG", 204, (204+(int(RT_test)-2)-(length_of_edit-1))) #determines if PAM is within 
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
                protosp_pe_b, PAM_site_b, guide_1_b, cas9_hp, edit_finish_b, pbs_finish_b, u6_term, pegrna_bot_finish, pe3_cut_distance_b, pe3_pe_finish_b, pe3_guide_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b, pbs_bot_final, DNA_wt_print, DNA_edit_print = pegrna_bot(DNA, Edit_site, Edit_nucl, edit_len_b, pbs_len_b, cas9_hp, u6_term, guide_api_request_out, input_position, design_alg)
                
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
            row_b = [pegRNA_name, chrom_num, input_position, editing_product, strand_b, pegrna_bot_finish, protosp_pe_b, edit_finish_b, edit_len_b, pbs_finish_b, pbs_bot_final, pe3_pe_finish_b, pe3_cut_distance_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b, seq_f_finish, seq_r_finish, notes, DNA_wt_print, DNA_edit_print]
            writer.writerow(row_b) #writes new row eachtime    
    
    if organism == "hg38":
        analysis_button = Button(root, text='Analyze', command=analyze_file).grid(row=10, column=3, columnspan=1)
    elif organism == "sacCer3":
        analysis_button = Button(root, text='Analyze', command=analyze_file).grid(row=10, column=3, columnspan=1)
        
    csvfile.close()
    return output_file_name

"""Circos Image Generation"""
#def open_circos():
#    global output_file_name
#    global circos_image #dont seem to need to declare image before global
#    top = Toplevel()
#    top.title('Analysis of ' +str(output_file_name))
#    lbl = Label(top, text='Hello?').pack()
#    circos_image = ImageTk.PhotoImage(Image.open("Low_quality_example.png"))
#    circos_label = Label(top, image=circos_image).pack()
#    close_button = Button(top, text='close', command=top.destroy).pack()


"""This is the Circos Analyze Functions"""
def analyze_file():
    global circos_image #dont seem to need to declare image before global
    
    def polar2xy(r, theta):
        return np.array([r*np.cos(theta), r*np.sin(theta)])
    
    def IdeogramArc(start=0, end=60, radius=1.0, width=0.2, ax=None, color=1):#color=(1,0,0) changing color here doesnt seem to affect
        # start, end should be in [0, 360]
        if start > end:
            start, end = end, start
        start *= np.pi/180.
        end *= np.pi/180.
        opt = 4./3. * np.tan((end-start)/ 4.) * radius
        inner = radius*(1-width)
        verts = [
            polar2xy(radius, start),
            polar2xy(radius, start) + polar2xy(opt, start+0.5*np.pi),
            polar2xy(radius, end) + polar2xy(opt, end-0.5*np.pi),
            polar2xy(radius, end),
            polar2xy(inner, end),
            polar2xy(inner, end) + polar2xy(opt*(1-width), end-0.5*np.pi),
            polar2xy(inner, start) + polar2xy(opt*(1-width), start+0.5*np.pi),
            polar2xy(inner, start),
            polar2xy(radius, start),
            ]
    
        codes = [Path.MOVETO,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.LINETO,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CLOSEPOLY,
                 ]
    
        if ax == None:
            return verts, codes
        else: #color infor for circle
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor=color, edgecolor='white', lw=2)
            ax.add_patch(patch) #remove "+(value)" to remove color distortion, +(0.5,), +(0.4,)
    
    """This is for chords between seperate sections""" # should not need to modify to much, just determines chord draw.
    def ChordArc(start1=0, end1=60, start2=180, end2=240, radius=1.0, chordwidth=0.7, ax=None, color=1):#color=(1,0,0)
        # start, end should be in [0, 360)
        if start1 > end1:
            start1, end1 = end1, start1
        if start2 > end2:
            start2, end2 = end2, start2
        start1 *= np.pi/180.
        end1 *= np.pi/180.
        start2 *= np.pi/180.
        end2 *= np.pi/180.
        opt1 = 4./3. * np.tan((end1-start1)/ 4.) * radius
        opt2 = 4./3. * np.tan((end2-start2)/ 4.) * radius
        rchord = radius * (1-chordwidth)
        verts = [
            polar2xy(radius, start1),
            polar2xy(radius, start1) + polar2xy(opt1, start1+0.5*np.pi),
            polar2xy(radius, end1) + polar2xy(opt1, end1-0.5*np.pi),
            polar2xy(radius, end1),
            polar2xy(rchord, end1),
            polar2xy(rchord, start2),
            polar2xy(radius, start2),
            polar2xy(radius, start2) + polar2xy(opt2, start2+0.5*np.pi),
            polar2xy(radius, end2) + polar2xy(opt2, end2-0.5*np.pi),
            polar2xy(radius, end2),
            polar2xy(rchord, end2),
            polar2xy(rchord, start1),
            polar2xy(radius, start1),
            ]
    
        codes = [Path.MOVETO,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 ]
    
        if ax == None:
            return verts, codes
        else: #color info for chords
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor=color+(0.5,), edgecolor=color+(0.5,), lw=0.3)
            ax.add_patch(patch) #modify  "+(value)" to manipulate color, +(0.5,), +(0.4,)
    
    """Plots chord diagram"""
    def chordDiagram(ax, colors=None, width=0.125, pad=0.0, chordwidth=0.7): #width is thickness of outer arc, pad is spacing between arcs, chordwith determines curvature of chords
        global valid_pegrna_sum
        
        if organism == 'sacCer3': #yeast genome dictionary
            chr_dict = {'chr1':230218,'chr2':813184,'chr3':316620,'chr4':1531933,
                              'chr5':576874,'chr6':270161,'chr7':1090940,'chr8':562643,
                              'chr9':439888,'chr10':745751,'chr11':666816,'chr12':1078177,
                              'chr13':924431,'chr14':784333,'chr15':1091291,'chr16':948066,'pegRNAs':600000}
        
        elif organism == 'hg38': #human genome dictionary
            chr_dict = {'chr1':248956422,'chr2':242193529,'chr3':198295559,'chr4':190214555,'chr5':181538259,
                          'chr6':170805979,'chr7':159345973,'chr8':145138636,'chr9':138394717,'chr10':133797422,
                          'chr11':135086622,'chr12':133275309,'chr13':114364328,'chr14':107043718,'chr15':101991189,
                          'chr16':90338345,'chr17':83257441,'chr18':80373285,'chr19':58617616,'chr20':64444167,
                          'chr21':46709983,'chr22':50818468,'chrx':156040895,'chry':57227415,'pegRNAs':145138636}
        
        chrom_names_list = list(chr_dict.keys())
        chrom_sizes_list = list(chr_dict.values())
        dictionary_sum = np.sum(chrom_sizes_list, dtype=np.int64).astype(float)
        pegRNA_degrees_start =  ((dictionary_sum-chr_dict['pegRNAs'])/dictionary_sum)*360
        pegRNA_degree_range = 360 - pegRNA_degrees_start
        #print(pegRNA_degree_range)
    
        ax.set_xlim(-1.1, 1.1) #sets axes scale, X
        ax.set_ylim(-1.1, 1.1) #sets axes scale, Y
    
        if colors is None: #color library 
            colors = [(0.40234375,0,0.05078125), (0.5977,0,0.0507), (0.79296875,0.09375,0.11328125), (0.93359375,0.23046875,0.171875),(0.98046875,0.4140625,0.2890625), (0.49609375,0.15234375,0.015625),
                      (0.6484375,0.2109375,0.01171875), (0.84765625,0.28125,0.00390625), (0.94140625,0.41015625,0.07421875),(0.98828125,0.55078125,0.234375), (0,0.265625,0.10546875), (0,0.42578125,0.171875),
                      (0.13671875,0.54296875,0.26953125), (0.25390625,0.66796875,0.36328125),(0.453125,0.765625,0.4609375), (0.03125,0.1875,0.41796875),(0.03125,0.31640625,0.609375), (0.12890625,0.44140625,0.70703125),
                      (0.2578125,0.5703125,0.7734375),(0.41796875,0.6796875,0.8359375), (0.24609375,0,0.48828125), (0.328125,0.15234375,0.55859375), (0.4140625,0.31640625,0.63671875), (0.5,0.48828125,0.7265625),
                      (0.03125,0.31640625,0.609375)]
            if len(chrom_sizes_list) > 25: #limit number of rows
                print('x is too large! Use x smaller than 10')
    
        #y determines proportion of arcs based off of chromosome sizes 
        y = chrom_sizes_list/np.sum(chrom_sizes_list, dtype=np.int64).astype(float) * (360 - pad*len(chrom_sizes_list)) #determines the proportion of the circle this takes up, 360 degrees- pad*number of arcs (rows)
        
        pos = {}
        arc = []
        nodePos = []    
        start = 0 #initial start position
    
        for i in range(len(chrom_sizes_list)): #x is sum of rows in data set, for loop runs and appends ARCS together
            end = start + y[i] #this determines arc size, use y[i]*0.5 to set to half circle
            arc.append((start, end)) #this looks like it adds new arcs to plot
            angle = 0.5*(start+end) #rotation of text
            
            if -30 <= angle <= 210:
                angle -= 90 #determines angle of text labels, 90 is tangent to circle, for top left half of circle
            else:
                angle -= 270 #determines angle of text labels, 90 is tangent to circle, for bottom right half of circle
            
            #######
            nodePos.append(tuple(polar2xy(0.93, 0.5*(start+end)*np.pi/180.)) + (angle,)) #controls distance out labels are places label for each arc
            #######
            
            start = end + pad #pad is distance between seperate arcs, start of new arc is defined here it appears
    
        for i in range(len(chrom_sizes_list)): #x is sum of rows in data set
            start, end = arc[i] #this detmines the arc position, essential
            #######
            IdeogramArc(start=start, end=end, radius=1.0, ax=ax, color=colors[i], width=width) #controls arc position, modifier can be added to end seq for chr. length
            #######
        """Analysis of excell output"""
        exe_path = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1]) #is able to return correct path to executable
        output_file_name_path = exe_path.replace('','') + '/' +output_file_name
        with open(output_file_name_path) as csvfile: #opens pinecone output
            readCSV = csv.reader(csvfile, delimiter=',')
            data1 = [list(row) for row in readCSV]
            #print(data)
            parameters1 = data1
            start2 = 340
            end2 = start2+0.5
            valid_pegrna_sum = 0
            
            for csv_add_value in range(1,len(parameters1)):
                pegRNA_status = str(parameters1[csv_add_value][5])
                if pegRNA_status == '':
                    none=0
                elif pegRNA_status == 'No Valid pegRNA':
                    none=0
                else:
                    valid_pegrna_sum += 1
            plotted_pegrna_count = 0
            print("Plotting " + str(valid_pegrna_sum) + " pegRNAs")
            for csv_index in range(1,len(parameters1)):
                pegRNA_name = str(parameters1[csv_index][0])
                pegRNA_chr = str(parameters1[csv_index][1]).lower()
                pegRNA_position = int(parameters1[csv_index][2])
                edit_result = str(parameters1[csv_index][3]).upper()
                pegRNA_status = str(parameters1[csv_index][5])
                
                if pegRNA_status == '':
                    none=0
                elif pegRNA_status == 'No Valid pegRNA':
                    none=0
                else:
                    input_key = "chr" + str(pegRNA_chr)
                    selection = chrom_names_list.index(input_key)
                    sum_total_chr = 0
                    sizes = 0
                    key_selection = chrom_names_list[0]
                    plotted_pegrna_count += 1
                    
                    q=0 #placeholder for q
                    for q in range(0,selection): #this is noninclusive of selection
                        key_selection = chrom_names_list[q] #seems to mess up here
                        sizes = chr_dict[key_selection]
                        sum_total_chr += sizes
                       
                    sum_total_chr_edit = sum_total_chr + pegRNA_position #i think this calc should happen outside this loop
                    edit_degrees = (sum_total_chr_edit/dictionary_sum)*(360-pad*(q)) #need modify pad term
                
                    width1 = 0.5 #width of chords connecting to chromosome
                    start1 = edit_degrees-width1*0.5
                    end1 = edit_degrees+width1*0.5
        
                    if pegRNA_degree_range/valid_pegrna_sum > 0.5: #this controls sizes of chords into 'pegRNA' node for small data sets              
                        width2 = 0.5
                        centering_term = pegRNA_degrees_start + (pegRNA_degree_range-(valid_pegrna_sum*0.5))*0.5
                        start2 = centering_term+plotted_pegrna_count*width2 #makes steps proportion to width
                        end2 = start2+width2
                        #print('A')
                    else: #this controls sizes of chords into 'pegRNA' node 
                        width2 = (pegRNA_degree_range/valid_pegrna_sum) #width will get smaller with large data sets
                        #print('B')
                        start2 = pegRNA_degrees_start+plotted_pegrna_count*width2 #makes steps proportion to width
                        end2 = start2+width2
                    
                    label_radius = 1.01
                    if 0 <= edit_degrees <= 90:
                        ax.annotate(pegRNA_name, polar2xy(label_radius, (edit_degrees)*np.pi/180.),horizontalalignment='left', verticalalignment='bottom', rotation=edit_degrees) #pegRNA name, radius, position, angle i think
                    elif 90 < edit_degrees <=180:
                        ax.annotate(pegRNA_name, polar2xy(label_radius, (edit_degrees)*np.pi/180.),horizontalalignment='right', verticalalignment='bottom', rotation=edit_degrees-180)
                    elif 180 < edit_degrees <=270:
                        ax.annotate(pegRNA_name, polar2xy(label_radius, (edit_degrees)*np.pi/180.),horizontalalignment='right', verticalalignment='top', rotation=edit_degrees-180)
                    elif 270 < edit_degrees <=359: #315
                        ax.annotate(pegRNA_name, polar2xy(label_radius, (edit_degrees)*np.pi/180.),horizontalalignment='left', verticalalignment='top', rotation=edit_degrees)
                    else:
                        ax.annotate(pegRNA_name, polar2xy(label_radius, (edit_degrees)*np.pi/180.),horizontalalignment='right', verticalalignment='center', rotation=edit_degrees-180) #pegRNA name, radius, position, angle i think
                    
                    if "D" in edit_result:
                        chord_color = (0.79296875,0.09375,0.11328125) #Red
                    elif "I" in edit_result:
                        chord_color = (0,0.42578125,0.171875) #Green
                    else:
                        chord_color = (0.03125,0.31640625,0.609375) #Blue,  (0.03125,0.31640625,0.609375)   
                    
                    #############
                    ChordArc(start1, end1, start2, end2, radius=1.-width, color=chord_color, chordwidth=chordwidth, ax=ax)
                    #############
                    
        #print(nodePos)
        return nodePos 
    
    """Circos Figure Plotting"""
    fig = plt.figure(figsize=(8,8)) #plotsize 10x10 works well
    ax = plt.axes([0,0,1,1])
    nodePos = chordDiagram(ax)
    ax.axis('off') #axes around plot
    prop = dict(fontsize=16*0.8, ha='center', va='center')
    if organism == "sacCer3":
        nodes_chr = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X','XI','XII','XIII','XIV','XV','XVI','pegRNAs']
    elif organism == "hg38":
        nodes_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','pegRNAs'] # Names of labels
    
    for i in range(len(nodes_chr)): # number here is number of labels, 
        ax.text(nodePos[i][0], nodePos[i][1], nodes_chr[i], rotation=nodePos[i][2], **prop, color='white')
    
    
    #Image Names
    exe_path = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1]) #is able to return correct path to executable
    output_file_name_path = exe_path.replace('','') + '/'
    high_res_display = output_file_name_path + "High_res_" + output_file_name_no_ext + ".svg"
    low_res_display = output_file_name_path + "Low_res_" + output_file_name_no_ext + ".png"
    
    #Save Images
    plt.savefig(high_res_display, format="svg",transparent=True, bbox_inches='tight', pad_inches=0.02) #1000dpi looks good
    plt.savefig(low_res_display, dpi=75, transparent=True, bbox_inches='tight', pad_inches=0.02) #1000dpi looks good
    
    
    """tkinter function to display plot"""
    top = Toplevel()
    top.title('Analysis of ' +str(output_file_name_no_ext))
    lbl = Label(top, text=str(output_file_name_no_ext)+" has "+str(valid_pegrna_sum)+" valid pegRNAs").pack()
    circos_image = ImageTk.PhotoImage(Image.open(low_res_display))
    circos_label = Label(top, image=circos_image).pack()
    close_button = Button(top, text='close', command=top.destroy).pack()

"""UI Startup Buttons and Images"""

ui_title = Label(root, text="PINE CONE: Creator Of New Edits", font=('helvetica', 18, 'bold')).grid(row=0, column=0, columnspan=4, pady=2)
ui_1 = Label(root, text="1: Select Organism", font=('helvetica', 12, 'bold')).grid(row=2, column=0, columnspan=1, pady=2)
ui_2 = Label(root, text="2: Select Design", font=('helvetica', 12, 'bold')).grid(row=2, column=1, columnspan=2, pady=2)
ui_3 = Label(root, text="3: Select Input Files", font=('helvetica', 12, 'bold')).grid(row=6, column=0, columnspan=1, pady=2)
ui_4 = Label(root, text="4: Output File", font=('helvetica', 12, 'bold')).grid(row=6, column=1, columnspan=2, stick=W, pady=2)
ui_5 = Label(root, text="Output Name:").grid(row=7, column=1, columnspan=2, stick=W, pady=2)

place_holder_csv_selected = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=7, column=0, columnspan=3, sticky=W)
place_holder_dna_input = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=8, column=0)
place_holder_dna_selected = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=10, column=1, columnspan=3)
place_holder_output_status = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=7, column=1, columnspan=3)

"""Organism Dropdown menu list of organisms, if oganis list updated, also update comboclick function"""

organism_options = ["Select Organism",
                    "Manual (.txt)",
                    "Plasmid (.txt)",
                    "Human (hg38)", 
                    "Yeast (S288C)",
                    "Mouse (mm10)",
                    "Zebrafish (danRer7)",
                    "Roundworm (ce11)",
                    "Fruitfly (dm6)"]




"""Guide preference Options"""

guide_option = ["Select Protospacer Preference",
                    "Nearest Protospacers to Edit",
                    "High Specificity (Requires Ref.Genome)"]

"""User Interface"""
organism_clicked = StringVar()
organism_clicked.set(organism_options[0])
#starting image for organism logo
starting_query_logo_path = os.path.join(logos_path,'Query_logo.png')
image1 = ImageTk.PhotoImage(Image.open(starting_query_logo_path))
label2 = Label(image=image1)
label2.image = image1
label2.grid(row=4, column=0, rowspan=2)

"""Dropdown Box for Organism Selection"""
myCombo = ttk.Combobox(root, value=organism_options) #dropdown menu
myCombo.current(0)
myCombo.bind("<<ComboboxSelected>>", comboclick) #binds combo, (thing, action), comboclick is function
myCombo.grid(row=3, column=0, padx=10)

"""Dropdown Box for Guide Design"""
myCombo1 = ttk.Combobox(root, value=guide_option, width=25) #dropdown menu
myCombo1.current(0)
myCombo1.bind("<<ComboboxSelected>>", comboclick1) #binds combo, (thing, action)
myCombo1.grid(row=3, column=1, columnspan=3)


"""Select File Button"""
my_btn = Button(root, text="Select Edit File (.csv)", command=select_file).grid(row=7, column=0, pady=2)

"""Banner Imadge"""
center_logo_path = os.path.join(logos_path,"Center_logo.png")
centerimage = ImageTk.PhotoImage(Image.open(center_logo_path))
centerlabel = Label(image=centerimage)
centerlabel.grid(row=1,column=0,columnspan=5, pady=5)

"""Labeling of buttons fg = foregroun and bg = background, can use Hex color code"""    
label_instruction = Label(root, text="Once file(s) are selected click 'Run PINE CONE'").grid(row=9, column=1, columnspan=3, sticky=W)
mybutton = Button(root, text="Run PINE CONE", command=run_button, fg='red', font=('helvetica', 12, 'bold')).grid(row=10, column=1, columnspan=3)

"""Prime Editing Strategy Editing Buttons"""
for index, (text, mode) in enumerate(MODES):
    editing_strategy_radiobutton = Radiobutton(root, text=text,variable=strategy_input, value=mode, command=lambda: clicked(strategy_input.get()))
    editing_strategy_radiobutton.grid(row=4,column=(1+index))

"""Output File Entry"""
output_file_entry = Entry(root)
output_file_entry.grid(row=7, column=2, columnspan=2, stick=W, pady=2)

"""Editing Strategy Image"""
def clicked(value):
    global image_label_strategy
    global strategy_selected
    image_label_strategy.grid_forget()
    strategy_selected = value
    strategy_icons = {'PE2_a':'PE2_Logo.png', 'PE3_a':'PE3_Logo.png', 'PE3B_a':'PE3B_Logo.png'}    
    strategy_logo_path = os.path.join(logos_path,strategy_icons[value])
    image_strategy = ImageTk.PhotoImage(Image.open(strategy_logo_path))
    image_label_strategy = Label(image=image_strategy)
    image_label_strategy.image = image_strategy
    image_label_strategy.grid(row=5, column=1, columnspan=3)
    print(strategy_selected)

root.mainloop()
print("Bye")