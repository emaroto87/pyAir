# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 15:27:25 2025

@author: U69432
"""

## Menu tipo windows ---------------------------------------------------------- <== MENU TIPO WINDOWS
import tkinter as tk
from tkinter.messagebox import showinfo
from ttkbootstrap import Window

from common.filedialog import askOpenOP2
from common.filedialog import askOpenGroupFile
from common.filedialog import askSaveAsExcel


def gui_windows():
    # Main Menu window
    # w = tk.Tk()
    # w.title('Prueba Menu')
    # w.geometry("600x300")
    w = Window(
        themename='flatly',
        # size = (600,600),
        title = 'NastranOP2Reader v1.0.2',
        resizable = None,
        scaling = False,
        )
    
    
    # Frame para colocar botones y opciones

    
    def show_info(wdg):
        showinfo(
            title = 'Result',
            message = wdg.get()
            )
    
    frame_load = tk.Frame(w)
    frame_load.pack(pady=20, padx=20)
    
    frame_rdbt = tk.Frame(w)
    frame_rdbt.pack(pady=50)

    
    # Boton de cargar OP2
    def select_op2_file():
        path = askOpenOP2()
        op2_path.delete(0,tk.END)
        op2_path.insert(0,path) 

    op2_lb = tk.Label(
        master = frame_load,
        text="OP2 file",
        width = 10,
        anchor="w"           
        )
    op2_lb.grid(row=0,column=0,padx=5,pady=5, sticky='w')
   
    op2_path = tk.Entry(
        master = frame_load,
        width=50
        )
    op2_path.grid(row=0,column=1, sticky='w')
    # Nota: si encadenas esta linea con la anterior, te dará el error de que
    # no encuentra la función delete() dentro de "select_op2_file". Esto es 
    # porque si lo metes entero, no guardas el widget, guardas el resultado de
    # colocarlo con la función grid, que es un None.
    
    op2_btn = tk.Button(
        master = frame_load,
        text = '\U0001F4C2',
        command = select_op2_file
        )
    op2_btn.grid(row=0,column=2, sticky='w')
    
        
    def select_xls_file():
        path = askSaveAsExcel()
        xls_path.delete(0,tk.END)
        xls_path.insert(0,path)
    
    xls_lb = tk.Label(
        master = frame_load,
        text="Excel file", 
        width = 10, 
        anchor="w"
        )
    xls_lb.grid(row=1,column=0,padx=5,pady=5, sticky='w')
    
    xls_path = tk.Entry(
        master = frame_load,
        width=50
        )
    xls_path.grid(row=1,column=1, sticky='w')

    xls_btn = tk.Button(
        master = frame_load,
        text = '\U0001F4C2',
        command = select_xls_file
        )
    xls_btn.grid(row=1,column=2, sticky='w')
    
    def select_set_file():
        path = askOpenGroupFile()
        set_path.delete(0,tk.END)
        set_path.insert(0,path)
    
    set_lb = tk.Label(
        master = frame_load,
        text="Set file", 
        width = 10, 
        anchor="w"
        )
    set_lb.grid(row=2,column=0,padx=5,pady=5, sticky='w')
    
    set_path = tk.Entry(
        master = frame_load,
        width=50
        )
    set_path.grid(row=2,column=1, sticky='w')

    set_btn = tk.Button(
        master = frame_load,
        text = '\U0001F4C2',
        command = select_set_file
        )
    set_btn.grid(row=2,column=2, sticky='w')
    
    
    # OPTIONS:
        
    frame_rdbt = tk.Frame(
        master = frame_load
        )
    frame_rdbt.grid(pady=50)
    
    # RadioButton : results options
    opt1 = tk.StringVar(value="")
    frame_results = tk.LabelFrame(
        master = frame_rdbt,
        text = 'Result Options'
        )
    frame_results.grid(row=0,column=0,padx=20,sticky='nw')
    for k,v  in __Result_Labels_Links.items():
        tk.Radiobutton(
            master = frame_results,
            text = k,
            variable = opt1,
            value = k,
            anchor = 'w',
            justify ='left'
            ).grid(sticky='w')
    
    # RadioButton : envelope_options
    opt2 = tk.StringVar(value="")

    frame_env = tk.LabelFrame(
        master = frame_rdbt,
        text = 'Envelope Options'
        )
    frame_env.grid(row=0,column=1,padx=20,sticky='nw')
    
    for env_opt in ENVELOPE_OPTIONS:
        tk.Radiobutton(
            master = frame_env,
            text = env_opt,
            variable = opt2,
            value = env_opt,
            anchor = 'w',
            justify ='left'
            ).grid(sticky='w')
    

    # RadioButton : Generate Summary
    opt3 = tk.StringVar(value="")
     
    frame_summary = tk.LabelFrame(
         master = frame_rdbt,
         text = 'Group Summary'
         )
    frame_summary.grid(row=0,column=3,padx=20,sticky='nw')
     
    yes_rbtn = tk.Radiobutton(
                 master = frame_summary ,
                 text = 'Yes',
                 variable = opt3,
                 value = 'Yes',
                 anchor = 'w',
                 justify ='left'
                 )
    yes_rbtn.grid(sticky='w')
    no_rbtn = tk.Radiobutton(
                 master = frame_summary ,
                 text = 'No',
                 variable = opt3,
                 value = 'No',
                 anchor = 'w',
                 justify ='left'
                 )
    no_rbtn.grid(sticky='w')
    
    w.mainloop()
## Menu tipo windows ---------------------------------------------------------  <== MENU TIPO WINDOWS