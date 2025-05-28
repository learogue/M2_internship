#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------------------------
# R script : Main script to call functions
# Auteur  : Léa ROGUE
# Date    : 31-03-2025
# Description : This python script show graphical interface to choose options to process and integrate scRNAseq.
# ----------------------------------------------------------------------------------------------------------------------------------------

# Load packages
import os
import re
import tkinter as tk
from tkinter import messagebox
import tkfilebrowser
from create_anndata_object import *
from apply_filters import *
from merge import *

# User chooses the processing directory
root = tk.Tk()
root.withdraw()  # To not display the main window
messagebox.showinfo('Information', 'Choose one directory for processing data or create one (icon in top right, write the name and enter). Click on OK to continue.')
processing_dir = tkfilebrowser.askopendirname(initialdir = '../../results/sc_results', title = 'Choose working directory or create it')
root.destroy()

# Create the user-specified subdirectories for Objects and Plots
if not os.path.isdir(processing_dir + '/Objects'):
    os.mkdir(processing_dir + '/Objects')
if not os.path.isdir(processing_dir + '/Plots'):
    os.mkdir(processing_dir + '/Plots')

# Create the graphical interface
root = tk.Tk()
root.geometry('300x280')
root.title('Select')

# Variables for the checkboxes
var1 = tk.BooleanVar()
var2 = tk.BooleanVar()
var3 = tk.BooleanVar()
var4 = tk.BooleanVar()
var5 = tk.BooleanVar()
var6 = tk.BooleanVar()

# Function called when the button is clicked
selected_options = []
def on_button_click():
    if var1.get():
        selected_options.append(1)
    if var2.get():
        selected_options.append(2)
    if var3.get():
        selected_options.append(3)
    if var4.get():
        selected_options.append(4)
    if var5.get():
        selected_options.append(5)
    if var6.get():
        selected_options.append(6)
    root.destroy()  # Close the main window

# Create the checkboxes and labels
label = tk.Label(root, text = 'Select what you want to do:')
label2 = tk.Label(root, text = 'Creating objects:')
label3 = tk.Label(root, text ='Manipulate objects')
label4 = tk.Label(root, text = 'Analyze objects')
checkbox2 = tk.Checkbutton(root, text = '10X format', variable = var2)
checkbox3 = tk.Checkbutton(root, text = 'Apply filters', variable = var3)
checkbox4 = tk.Checkbutton(root, text = 'Merge objects', variable = var4)
checkbox5 = tk.Checkbutton(root, text = 'Create umaps', variable = var5)
checkbox6 = tk.Checkbutton(root, text = 'Find differentially expressed genes per cluster', variable = var6)

# Position the checkboxes and labels
label.pack()                # Add the label to the window
label2.pack(ancho = tk.W)    # Use anchor to align the button
checkbox2.pack()
label3.pack(anchor = tk.W)
checkbox3.pack()
checkbox4.pack()
label4.pack(anchor = tk.W)
checkbox5.pack()
checkbox6.pack()

# Button to get the selected options
btn_get_selected = tk.Button(root, text = 'Select', command = on_button_click)
btn_get_selected.pack()

# Start the main loop
root.mainloop()

# Processing for 10X
if 2 in selected_options:
    # Part 1: Create Anndata objects for 10X
    root = tk.Tk()
    root.withdraw()
    messagebox.showinfo('Information', 'Choose one or more 10X folder(s) to create an object')
    l_dir = list(tkfilebrowser.askopendirnames(initialdir = '../../data/scrnaseq_data/', title = 'Choose directories for 10X'))

    # Create Anndata objects for each specified 10X directory
    for dir in l_dir:
        create_anndata_object(re.split(r'[/\\]', dir)[-1], '10X', processing_dir)
    messagebox.showinfo('Information', 'Object(s) creation finished')
    root.destroy()

# Apply filters section
if 3 in selected_options:
    root = tk.Tk()
    root.withdraw()
    messagebox.showinfo('Information', 'Choose one or more object(s) to apply filters')
    l_obj = list(tkfilebrowser.askopenfilenames(initialdir = processing_dir+'/Objects/Objects_ori', title = 'Choose objects to apply filters'))

    # Apply filters to each specified Anndata object
    for obj in l_obj:
        apply_filters(re.split(r'[/\\]', obj)[-1], processing_dir)
    messagebox.showinfo('Information', 'Finished')
    root.destroy()

# Merge section
if 4 in selected_options:
    root = tk.Tk()
    root.withdraw()
    messagebox.showinfo('Information', 'Step of merging some object')
    l_path_obj = list(tkfilebrowser.askopenfilenames(initialdir = processing_dir+'/Objects/Objects_filtered/', title = 'Choose objects for merging'))

    # Take only the object name
    l_obj = []
    for obj in l_path_obj:
        l_obj.append(re.split(r'[/\\]', obj)[-1])

    # Merge the specified objects
    merge(l_obj, processing_dir)
    messagebox.showinfo('Information', 'Merging finished')
    root.destroy()

