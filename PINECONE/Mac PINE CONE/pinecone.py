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
import numpy as np #used in tm caclulation
#########

print("Hello")
#When running as source app_path = 'pinecone/Contents/MacOS'. when in App app_path = 'pinecone.app/Contents/MacOS'
root = Tk()
root.title('PINE CONE')
exe_path = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1]) #is able to return correct path to executable
logos_path = os.path.join(exe_path,'logos')
app_path = 'pinecone.app/Contents/MacOS' #pinecone.app/Contents/MacOS
program_logo_path = os.path.join(logos_path,'pineconelogo.icns')
root.iconbitmap(program_logo_path) #logo
root.geometry("510x550") #Size of UI window
#frame = LabelFrame(root, text="Frame", padx=5, pady=5, bg='#f0f0f0')
#frame.grid(row=0, column=0)

"""OPENS CSV INPUT PARAMETERS"""
input_file_path = "aaaaaa"
input_file_path_text = "zzzzz"
parameters = "bbbbbb"
organism = "null"
analysis_button = "null"
strategy_selected = "PE2_a"

"""UI Stuff"""
"""Organism Dropdown menu, selection dictionary defines organism searched"""
def comboclick(event):
    global label2
    global organism
#    myLabel1=Label(root, text=myCombo.get())
#    myLabel1.grid(row=0, column=0)
    
    print(myCombo.get())
    
    label2.grid_forget()
    #dictionary converts text into general term, t be used elsewhere (for more complicated programs)
    selction_dictionary = {'Select Organism':'null',
                           'Manual (.txt)':'text',
                           "Plasmid (.txt)":'plasmid',
                           'Human (hg38)':'hg38', 
                           'Yeast (S288C)':'S288C',
                           'Mouse (mm10)':'mm10',
                           'Rat (rn6)':'rn6',
                           'Zebrafish (danRer11)':'danRer11',
                           'Roundworm (ce11)':'ce11',
                           'Fruitfly (dm6)':'dm6'}
    
    selection_intermediate = selction_dictionary[myCombo.get()] #just an intermediate to translate between dictionaries
    organism = selction_dictionary[myCombo.get()] #organism is name of reference genome as appears in URL
    #dictionary converts selection variable to image name
    organism_icons = {'null':'Query_logo.png', 
                      'hg38':'Human_logo.png', 
                      'S288C':'yeast_logo.png',
                      'mm10':'mouse_logo.png',
                      'text':'manual_text_logo.png',
                      'plasmid':'manual_text_logo.png',
                      'rn6':'rat_logo.png',
                      'danRer11':'zebrafish_logo.png',
                      'ce11':'roundword_logo.png',
                      'dm6':'fruitfly_logo.png'}
    
    image_to_display = os.path.join(logos_path, organism_icons[selection_intermediate])

    if myCombo.get() == "Manual (.txt)":
#        myLabel = Label(root, text = 'select text document').pack()
        my_btn = Button(root, text="Select DNA Sequence (.txt)", command=select_text_file).grid(row=8, column=0) #if manual selected creates file dialog for text input
        #print(myCombo.get())
    elif myCombo.get() == "Plasmid (.txt)":
        my_btn = Button(root, text="Select DNA Sequence (.txt)", command=select_text_file).grid(row=8, column=0)
    else:
        status_label_organism = "Organism: " + myCombo.get()
        #myLabel = Label(root, text = status_label_organism).grid(row=7, column=1, columnspan=3, sticky=W)

    image1 = ImageTk.PhotoImage(Image.open(image_to_display))
    label2 = Label(image=image1)
    label2.image = image1
    label2.grid(row=4, column=0)

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
image_label_strategy.grid(row=4, column=1, columnspan=3)



"""Select .CSV file"""
def select_file():
    global input_file_path
    root.filename = filedialog.askopenfilename(initialdir = "/", title="Select File", filetypes = (("csv files", "*.csv"),("all files","*.*")))
    csv_file_name = os.path.basename(root.filename)
    csv_file_display = "Edit File: " + csv_file_name
    my_label = Label(root, text=csv_file_display).grid(row=7, column=0, columnspan=1) #sticky=W
    print (root.filename)
    input_file_path = root.filename

"""If manual input selected, file select for .txt file"""
def select_text_file():
    global input_file_path_text
    root.filename_text = filedialog.askopenfilename(initialdir = "/", title="Select File", filetypes = (("text files", "*.txt"),("all files","*.*")))
    text_file_name = os.path.basename(root.filename_text)
    text_file_display = "DNA Sequence: " + text_file_name
    my_label = Label(root, text=text_file_display).grid(row=9, column=0, columnspan=1)#sticky=W
    input_file_path_text = root.filename_text
    #print(root.filename_text)

"""Main Function, Run button initiates this"""
def run_button():
    global parameters
    global organism
    global output_file_name
    global output_file_name_no_ext
    
    if output_file_entry.get() == "": #Sets default name of output if no output name is entered
        output_file_name_no_ext = "pinecone_output"
        output_file_name = output_file_name_no_ext + ".csv"
    else: #Sets output name to user input
        output_file_name_no_ext = output_file_entry.get()
        output_file_name = output_file_name_no_ext + ".csv"
        
    output_file_status = "'" + output_file_name + "' is complete"
    clicklabel = Label(root, text=output_file_status, fg='green', font=('helvetica', 12, 'bold'))
    clicklabel.grid(row=7, column=1, columnspan=3)
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
     
        
    def primer_tm_design(sequence): #enter dna sequence, generates primer with specific Tm, input seq should be 40bp or longer
        # format for primer requests, generates dictionary with primer information
        tm_set = 63  # Sets goal Tm
        tm = 0  # Starting Tm at begining of For-loop
        # dH is deltaH enthalpy dinucleotide dictionary
        dH = {'AA': -7.9, 'TT': -7.9, 'AT': -7.2, 'TA': -7.2,
              'CA': -8.5, 'TG': -8.5, 'GT': -8.4, 'AC': -8.4,
              'CT': -7.8, 'AG': -7.8, 'GA': -8.2, 'TC': -8.2,
              'CG': -10.6, 'GC': -9.8, 'GG': -8.0, 'CC': -8.0}
        # dS is entropy, -kcal/mol dinucleotide dictionary
        dS = {'AA': -22.6, 'TT': -22.6, 'AT': -20.4, 'TA': -21.3,
              'CA': -22.7, 'TG': -22.7, 'GT': -22.4, 'AC': -22.4,
              'CT': -21.6, 'AG': -21.6, 'GA': -22.2, 'TC': -22.2,
              'CG': -27.2, 'GC': -24.4, 'GG': -19.9, 'CC': -19.9}

        adjust_H = 0.1
        adjust_S = 2.8

        dH_sum = 0 + adjust_H
        dS_sum = 0 + adjust_S

        for base in range(0, len(sequence) - 2):
            dinucleotide = sequence[base:base + 2]
            dH_sum += dH[dinucleotide]
            dS_sum += dS[dinucleotide]
            R = 1.987  # universal gas constant
            M = 500 * 10 ** -9  # common molarity of primers
            tm_int = ((dH_sum) * 1000) / (dS_sum + R * np.log(M)) - 273.15
            tm = (tm_int + 16.6 * np.log(1)) - 15  # 15 is adjustment to roughly match values
            if tm >= tm_set:  # elongates primer until Tm is reached
                break

        primer_seq = sequence[0:base + 2]  # this is final primer sequence
        # print("Primer legnth: "+ str(len(primer_seq)))
        # print("Primer Seq: " + primer_seq)
        # print("Calculated TM is: " + str(tm))

        return primer_seq
    
    """DESIGNS POSITIVE STRAND PEGRNA"""
    def pegrna_top(DNA, Edit_site, Edit_nucl, edit_len_a, pbs_len_a, cas9_hp, u6_term):
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
            DNA_r = DNA[1+Edit_site:len(DNA)]
            DNA_edited = DNA_l + Edit_nucl + DNA_r

        length_of_edit = len(Edit_nucl) #computes length of Edit in function, used to modify PAM search and avoid distant PAMs
                
        dna_reverse = DNA[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
        pam_rev_search = dna_reverse.index("GG", 205, (205+(edit_len-2)-(length_of_edit-1))) #default 203,202+edit_len_a-(length_of_edit-1) #determines PAM site from reverse, to place as close to edit as possible , 203 (reverse) is left most PAM and 200+edit_len_a is edit range with 3bp buffer in RT 
        PAM_site = len(DNA)-pam_rev_search-2
        
        #PAM_site = DNA.index("GG", Edit_site-1, Edit_site+6)
        
        protosp_pe = DNA[PAM_site-21:PAM_site-1]
        guide_1 = "G" + DNA[(PAM_site-21):(PAM_site-1)]
    
        pbs_1 = DNA[(PAM_site-4-pbs_len_a):(PAM_site-4)]
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
        print("Strategy is: " +str(strategy_selected))
        
        if strategy_selected == 'PE2_a':
            pe3_cut_distance = "" #blank since PE2 selected
            pe3_pe_finish = "" #blank since PE2 selected
            pe3_guide = "" #blank since PE2 selected
            ### Oligos ###
            pe3_guide_top = "" #blank since PE2 selected
            pe3_guide_bot = "" #blank since PE2 selected
            
        elif strategy_selected == 'PE3_a':
            try: #attempts to find PE3 protospacer
                pe3_pam_site = DNA.index("CC", PAM_site+30,PAM_site+100) #default +40,+90#searches for PE3 guide PAM on opposite strand
                pe3_pe_start = DNA[int(pe3_pam_site)+3:int(pe3_pam_site)+23] #double check spacings
                pe3_cut_distance = int(pe3_pam_site)-int(PAM_site)+10 #
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_pe_start) 
                pe3_pe_finish = dna_rc #
                pe3_guide = "G" + pe3_pe_finish
                ### Oligos ###
                pe3_guide_top = "CACC" + pe3_guide
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_guide)
                pe3_guide_bot = "AAAC" + dna_rc
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
                pe3_pam_site_init = DNA_edited_reverse.index("CC", Edit_site+3,Edit_site+23) #default +40,+90#searches for PE3 guide PAM on opposite strand
                pe3_pam_site = len(DNA)-pe3_pam_site_init-2
                pe3_pe_start = DNA_edited[int(pe3_pam_site)+3:int(pe3_pam_site)+23] #double check spacings
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
        
        return protosp_pe, PAM_site, guide_1, cas9_hp, edit_finish, pbs_finish, u6_term, pegrna_top_finish, pe3_cut_distance, pe3_pe_finish, pe3_guide, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot     
    
    """DESIGNS BOTTOM STRAND pegRNA"""
    def pegrna_bot(DNA, Edit_site, Edit_nucl, edit_len_a, pbs_len_a, cas9_hp, u6_term):
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
            DNA_r = DNA[1+Edit_site:len(DNA)]
            DNA_edited = DNA_l + Edit_nucl + DNA_r
        
        dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA)
        dna_rcs = str(dna_rc)
        dna_reverse = dna_rcs[::-1]

        length_of_edit = len(Edit_nucl) #computes length of Edit in function, used to modify PAM search and avoid distant PAMs
        pam_rev_search = dna_reverse.index("GG", 204, (204+(edit_len-2)-(length_of_edit-1))) #determines if PAM is within 
        PAM_site_b = len(DNA)-pam_rev_search-2
        protosp_pe_b = dna_rcs[PAM_site_b-21:PAM_site_b-1]
        guide_1_b = "G" + dna_rcs[(PAM_site_b-21):(PAM_site_b-1)]
        pbs_1 = dna_rcs[(PAM_site_b-4-pbs_len_a):(PAM_site_b-4)] #defining PBS
        dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pbs_1)
        pbs_finish_b = dna_rc
        
        dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA_edited) #defining RT edit
        edit_1 = dna_rc[(PAM_site_b-4):(PAM_site_b-4+edit_len_a)]
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
                pe3_pam_site_b = dna_rcs.index("CC", PAM_site_b+30,PAM_site_b+100) #default +40,+90#searches for PE3 guide PAM on opposite strand
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
                
                pe3_pam_site_init = DNA_edited_reverse.index("CC", Edit_site+3,Edit_site+23) #PE3B defined based off of edit
                
                pe3_pam_site_b = len(DNA)-pe3_pam_site_init-2
                pe3_pe_start = DNA_edited[int(pe3_pam_site_b)+3:int(pe3_pam_site_b)+23] #double check spacings
                pe3_cut_distance_b = int(pe3_pam_site_b)-int(PAM_site_b)+10 #
                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(pe3_pe_start) 
                pe3_pe_finish_b = dna_rc #
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
    
        return protosp_pe_b, PAM_site_b, guide_1_b, cas9_hp, edit_finish_b, pbs_finish_b, u6_term, pegrna_bot_finish, pe3_cut_distance_b, pe3_pe_finish_b, pe3_guide_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b  
    
    """LOOP FOR GENERATING PEGRNA CSV FILE"""
    output_file_name_path = exe_path.replace(app_path,'') + output_file_name #sets write file path to same directory as app
    #output_file_name_path = output_file_name
    #.replace('Pinecone_test','output_dir/')

    with open(output_file_name_path,'w') as f1: #Creates CSV file
        writer=csv.writer(f1, delimiter=',',lineterminator='\n',) #assigns writer parameters, \t
        header = ['Name','Chromosome','Position', 'Edit', 'Strand', 'pegRNA','Protospacer' ,'RT Template','RT Length','Primer Sequence','PBS Length', 'PE3 or 3B Protospacer', 'PE3 or 3B Distance', 'peg_Guide_top', 'peg_Guide_bot', 'peg_Edit_top', 'peg_Edit_bot', 'PE3_Guide_top', 'PE3_Guide_bot', 'Seq_F', 'Seq_R', 'Notes'] #header of CSV file output
        writer.writerow(header) #creates header in file
        for index in range(1,len(parameters)):
            ###Input_CSV###
            pegRNA_name = str(parameters[index][0]) #Name of guide in input file
            chrom_num = str(parameters[index][1]) #chromosome of target locus
            input_position = int(parameters[index][2]) #specific Bp address of locus
            Edit_nucl = str(parameters[index][3]).upper() #Edit Specified in RT edit
            length_of_edit = len(parameters[index][3]) #used for modifying PAM site search calculation, avoid spurious ID of distant PAMs
            #Edit_site = int(209) #this is hard coded position of edit in seq returned from API
            edit_len = int(parameters[index][4]) #RT template length, also specifies editing range for PAM search
            edit_len_a = edit_len #explicitly keeps edit length as an integer
            pbs_len = int(parameters[index][5]) #Specified Primer length
            pbs_len_a = pbs_len #explicitly kepps PBS length as an integer
            notes = str(parameters[index][6]) #This is retained between input and output for note keeping
            #https://genome.ucsc.edu/goldenPath/help/api.html > get DNA seq from speficied chromosomes
#            position_start = str(input_position-210) #downstream of edit
#            position_end = str(input_position+210) #upstream of edit
#            url_req = str("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr" + chrom_num + ";start=" + position_start +";end=" + position_end) #URL of UCSC Genome browaer Hg38 API
            try:
                ##################### Organism selection #####################
                if organism == "null":
                    print("no organism selected")
                    DNA = "No Organim, No DNA"
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
                    #print(DNA)
                    
                elif organism == "plasmid": #need to add text wrap around
                    file_name = open(input_file_path_text, "r")
                    print(file_name.read)
                    print(file_name.readable()) #TRUE OR FALSE output
                    DNA_in = str(file_name.readline()).upper().replace(' ','').replace('\n','') #converts all to upper and removes spaces and new lines
                    Edit_site = int(209)
                    
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
                    #print(DNA)
                
                elif organism == "S288C": #Yeast Saccharomyces cerevisae S288C, YGD
                    position_start = str(input_position-209) #downstream of edit
                    position_end = str(input_position+210) #upstream of edit
                    url_req = str("https://www.yeastgenome.org/run_seqtools?&chr=" + chrom_num + "&start=" + position_start +"&end=" + position_end + "&rev=0")
                    api_request = requests.get(url_req) #Begins API request
                    api = json.loads(api_request.content)
                    DNA_in = api['residue']
                    DNA = DNA_in.upper().replace(' ','') #converts all to upper and removes spaces
                    Edit_site = int(209) #this is hard coded position of edit in seq returned from API
                    print(DNA) 
                else: #Human Homo sapiens Hg38 and other
                    position_start = str(input_position-210) #downstream of edit
                    position_end = str(input_position+210) #upstream of edit
                    url_req = str("https://api.genome.ucsc.edu/getData/sequence?genome=" + organism + ";chrom=chr" + chrom_num + ";start=" + position_start +";end=" + position_end) #URL of UCSC Genome browaer Hg38 API
                    api_request = requests.get(url_req) #Begins API request
                    api = json.loads(api_request.content)
                    DNA_in = api['dna']
                    DNA = DNA_in.upper().replace(' ','') #converts all to upper and removes spaces
                    Edit_site = int(209) #this is hard coded position of edit in seq returned from API
                    print(DNA)

                ##################### Test if pegRNA can be designed, Positive Strand #####################
                dna_reverse = DNA[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                pam_rev_search = dna_reverse.index("GG", 205, 205+((edit_len-2)-(length_of_edit-1))) #determines PAM site from reverse, to place as close to edit as possible 
                PAM_site = len(DNA)-pam_rev_search-2
    
                #PAM_site = DNA.index("GG", Edit_site-1, Edit_site+6) #OLD PAM SITE DEFINITION
                protosp_pe, PAM_site, guide_1, cas9_hp, edit_finish, pbs_finish, u6_term, pegrna_top_finish, pe3_cut_distance, pe3_pe_finish, pe3_guide, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot = pegrna_top(DNA, Edit_site, Edit_nucl, edit_len, pbs_len, cas9_hp, u6_term)
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
                    editing_product = str(DNA[209]) + "-to-" + str(Edit_nucl) #generates X-to-Y format for edit annotation
                
                print("PE3 Cut Distance: " + str(pe3_cut_distance))
                print("PE3 Protospacer: " + pe3_pe_finish)
                print("Forward Sequencing: " + seq_f_finish)
                print("Reverse Sequencing: " + seq_r_finish)
                print("Length of Edit " + str(index) + ": " + str(length_of_edit))
                strand_t = "Positive"
                #validation_top = "valid"
            except ValueError:
                pegrna_top_finish = "No Valid pegRNA"
                strand_t = "Positive"
                seq_f_finish = " "
                seq_r_finish = " "
                editing_product = " "
                protosp_pe = " "
                edit_finish = " "
                pbs_finish = " "
                pe3_pe_finish = " "
                pe3_cut_distance = " "
                peg_guide_top = " "
                peg_guide_bot = " "
                peg_edit_top = " "
                peg_edit_bot = " "
                pe3_guide_top = " "
                pe3_guide_bot = " "
                #validation_top = "error"
                print("No Valid positive strand pEGRNAs for " + str(index))
            row_t = [pegRNA_name, chrom_num, input_position, editing_product, strand_t, pegrna_top_finish, protosp_pe, edit_finish, edit_len, pbs_finish, pbs_len, pe3_pe_finish, pe3_cut_distance, peg_guide_top, peg_guide_bot, peg_edit_top, peg_edit_bot, pe3_guide_top, pe3_guide_bot, seq_f_finish, seq_r_finish, notes] #determines what values are writen to CSV file 
            writer.writerow(row_t) #writes new row eachtime
            #Negative strand pegRNAs below
            try:
                dna_rc, dna_c, dna_r, dna_a  = dna_rev_comp(DNA)
                dna_rcs = str(dna_rc)
                
                dna_reverse = dna_rcs[::-1] #reverses DNA for PAM search, PAM search done on reverse to favor guides close to edit
                pam_rev_search = dna_reverse.index("GG", 204, (204+(edit_len-2)-(length_of_edit-1))) #determines if PAM is within 
                PAM_site_b = len(DNA)-pam_rev_search-2
                protosp_pe_b, PAM_site_b, guide_1_b, cas9_hp, edit_finish_b, pbs_finish_b, u6_term, pegrna_bot_finish, pe3_cut_distance_b, pe3_pe_finish_b, pe3_guide_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b = pegrna_bot(DNA, Edit_site, Edit_nucl, edit_len, pbs_len, cas9_hp, u6_term)
                print("Bottom strand protospacer " + str(index) + ":" + protosp_pe_b)
                print("Bottom pEGRNA for " + str(index) + ":")
                print("pegRNA for "  + str(index) + ":")
                #pegrna_bot_finish = str(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term)
                print(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term)
                pegrna_bot_finish = str(guide_1_b + cas9_hp + edit_finish_b + pbs_finish_b + u6_term)

                ### Primers ###
                sequence = str(DNA[0:60])
                primer_seq = primer_tm_design(sequence)  # forward primer, tm set
                seq_f_finish = primer_seq

                dna_rc, dna_c, dna_r, dna_a = dna_rev_comp(DNA[len(DNA) - 60:len(DNA)])  # reverse primer, tm set
                sequence = str(dna_rc)
                print(sequence)
                primer_seq = primer_tm_design(sequence)
                seq_r_finish = primer_seq
                
                if "D" in Edit_nucl: #Determines if user is requesting deletion
                    editing_product = str(Edit_nucl)
                elif "I" in Edit_nucl:
                    editing_product = str(Edit_nucl)
                else:
                    editing_product = str(DNA[209]) + "-to-" + str(Edit_nucl) #generates X-to-Y format for edit annotation
                
                print("PE3 Cut Distance: " + str(pe3_cut_distance_b))
                print("PE3 Protospacer: " + pe3_pe_finish_b)
                print("Forward Sequencing: " + seq_f_finish)
                print("Reverse Sequencing: " + seq_r_finish)
                #validation_bot = "valid"
                strand_b = "Negative"
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
                pe3_pe_finish_b = " "
                pe3_cut_distance_b = " "
                peg_guide_top_b = " "
                peg_guide_bot_b = " "
                peg_edit_top_b = " "
                peg_edit_bot_b = " "
                pe3_guide_top_b = " "
                pe3_guide_bot_b = " "
                #validation_bot = "error"
                print("No Valid negative strand  pEGRNAs " + str(index))
            row_b = [pegRNA_name, chrom_num, input_position, editing_product, strand_b, pegrna_bot_finish, protosp_pe_b, edit_finish_b, edit_len, pbs_finish_b, pbs_len, pe3_pe_finish_b, pe3_cut_distance_b, peg_guide_top_b, peg_guide_bot_b, peg_edit_top_b, peg_edit_bot_b, pe3_guide_top_b, pe3_guide_bot_b, seq_f_finish, seq_r_finish, notes]
            writer.writerow(row_b) #writes new row eachtime

            if organism == "text":
                time.sleep(0)
            elif organism == "plasmid":
                time.sleep(0)
            else:
                time.sleep(0.5) #limits API request to 1 per half second, prevents request spamming
    

    csvfile.close()
    return output_file_name

"""UI Startup Buttons and Images"""

ui_title = Label(root, text="PINE CONE: Creator Of New Edits", font=('helvetica', 18, 'bold')).grid(row=0, column=0, columnspan=4, pady=2)
ui_1 = Label(root, text="1: Select Organism", font=('helvetica', 12, 'bold')).grid(row=2, column=0, columnspan=1, pady=2)
ui_2 = Label(root, text="2: Select Strategy", font=('helvetica', 12, 'bold')).grid(row=2, column=1, columnspan=2, pady=2)
ui_3 = Label(root, text="3: Select Input Files", font=('helvetica', 12, 'bold')).grid(row=5, column=0, columnspan=1, pady=2)
ui_4 = Label(root, text="4: Output File", font=('helvetica', 12, 'bold')).grid(row=5, column=1, columnspan=2, stick=W, pady=2)
ui_5 = Label(root, text="Output Name:").grid(row=6, column=1, columnspan=2, stick=W, pady=2)

place_holder_csv_selected = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=7, column=0, columnspan=3, sticky=W)
place_holder_dna_input = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=8, column=0)
place_holder_dna_selected = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=9, column=1, columnspan=3)
place_holder_output_status = Label(root, text=" ", font=('helvetica', 12, 'bold')).grid(row=7, column=1, columnspan=3)

"""Organism Dropdown menu list of organisms, if oganis list updated, also update comboclick function"""

organism_options = ["Select Organism",
                    "Manual (.txt)",
                    "Plasmid (.txt)",
                    "Human (hg38)", 
                    "Yeast (S288C)",
                    "Mouse (mm10)",
                    "Rat (rn6)",
                    "Zebrafish (danRer11)",
                    "Roundworm (ce11)",
                    "Fruitfly (dm6)"]


"""User Interface"""
organism_clicked = StringVar()
organism_clicked.set(organism_options[0])
#starting image for organism logo
starting_query_logo_path = os.path.join(logos_path,'Query_logo.png')
image1 = ImageTk.PhotoImage(Image.open(starting_query_logo_path))
label2 = Label(image=image1)
label2.image = image1
label2.grid(row=4, column=0)

"""PE Editing Strategy Radio Buttons"""


"""Dropdown Box for Organism Selection"""
myCombo = ttk.Combobox(root, value=organism_options)
myCombo.current(0)
myCombo.bind("<<ComboboxSelected>>", comboclick)#binds combo, (thing, action)
myCombo.grid(row=3, column=0)


"""Select File Button"""
my_btn = Button(root, text="Select Edit File (.csv)", command=select_file).grid(row=6, column=0, pady=2)

"""Banner Image"""
center_logo_path = os.path.join(logos_path,"Center_logo.png")
centerimage = ImageTk.PhotoImage(Image.open(center_logo_path))
centerlabel = Label(image=centerimage)
centerlabel.grid(row=1,column=0,columnspan=5, pady=5)

"""Labeling of buttons fg = foregroun and bg = background, can use Hex color code"""    
label_instruction = Label(root, text="Once file(s) are selected click 'Run PINE CONE'").grid(row=8, column=1, columnspan=3, sticky=W)
mybutton = Button(root, text="Run PINE CONE", command=run_button, fg='red', font=('helvetica', 12, 'bold')).grid(row=9, column=1, columnspan=3)
#fg='white', bg='red',
"""Prime Editing Strategy Editing Buttons"""
for index, (text, mode) in enumerate(MODES):
    editing_strategy_radiobutton = Radiobutton(root, text=text,variable=strategy_input, value=mode, command=lambda: clicked(strategy_input.get()))
    editing_strategy_radiobutton.grid(row=3,column=(1+index))

"""Output File Entry"""
output_file_entry = Entry(root)
output_file_entry.grid(row=6, column=2, columnspan=2, stick=W, pady=2)

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
    image_label_strategy.grid(row=4, column=1, columnspan=3)

root.mainloop()
print("Bye")