# -*- coding: utf-8 -*-
"""
Tkinter Button widget:
======================

Button(master, **options)
  text >> Texto para mostrar en el botón
  image >> Image para mostrar en el boton. Deberá de ser un objeto PhotoImage, BitmapImage o compatible
  height >> Alto del boton
  width >> Ancho del boton
  state >> Estado del boton --> [NORMAL, ACTIVE, DISABLED] Default NORMAL
  command >> Función a la que se asocia el click del boton
  background // bg >> Color de fondo del boton. String con color hexadecimal ("#000fff000) o color de sistema ('red')
  foreground // fg >> Color para texto del botón. String con color hexadecimal ("#000fff000) o color de sistema ('red')
  borderwidth // bd >> Ancho de borde
  relief >> FLAT, SUNKEN, RAISED, GROOVE, RIDGE. Necesario cambiar propiedad borderwidth
  anchor >> Controla donde se colocará el texto dentro del boton [N, NE, E, SE, S, SW, W, NW, CENTER] Default CENTER
  justify >> Controla como se justifica el texto. Default CENTER
  font >> Fuente para utilizar en el texto del botón --> "nombre tamaño weight style"


Para asociar varios botones con una sola funcion:
1. Crear una función que tome un argumento (numero entero)
2. Crear lista de botones y en command utilizar la función lambda

def callback(indx):
    pass # Do something

botones = []
for n in range(4):
    b = Button(root, text="Abrir", command=lambda idx=n : callback(idx))  # Necesario porque n cambia y no lo coge bien
    b.pack(side=LEFT)
    botones.append(b)

"""

from tkinter import Frame, Label, Tk, Button, Entry, StringVar, filedialog
from tkinter import LEFT, RIGHT, TOP, BOTTOM,  W, N, E, S, DISABLED, ACTIVE, YES, NO, X, Y, BOTH


# Form desing
class MyForm:
    def __init__(self, master):
        self.master = master
        self.master.title("Get channels")

        # Creamos 3 frames con botones y etiquetas
        texts = ["Digital Elevation Model:", "Flow Accumulation:", "Heads shapefile:", "Basins shapefile:"]
        self.entry_vars = [StringVar() for n in range(3)]
        for n in range(3):
            f = Frame(self.master)
            l = Label(f, width=19, anchor=W, justify=LEFT, text=texts[n])
            l.pack(side=LEFT, padx=3, expand=NO)
            tx = Entry(f, textvariable=self.entry_vars[n])
            tx.pack(side=LEFT, padx=3, expand=YES, fill=X)
            b = Button(f, text="Abrir", command=lambda idx=n: self.open_file(idx))
            b.pack(side=LEFT, padx=3, expand=NO)
            f.pack(side=TOP, padx=2, pady=2, expand=NO, fill=X, anchor=N)

        # Creamos otro frame con dos botones
        f = Frame(self.master)
        self.quit_btn = Button(f, width=10, text="Quit", command=self.master.quit).pack(side=RIGHT, padx=3, expand=NO)
        self.run_btn = Button(f, width=10, text="Run", command=self.run_program).pack(side=RIGHT, padx=3, expand=NO)
        f.pack(side=TOP, padx=2, pady=10, expand=NO, fill=X, anchor=N)

    def run_program(self):
        print("Archivos seleccionados:")
        for s_var in self.entry_vars:
            print(s_var.get())

    def open_file(self, idx):
        filename = filedialog.askopenfilename(filetypes=[('All files', '.pdf')])
        self.entry_vars[idx].set(filename)

root = Tk()
myapp = MyForm(root)
root.mainloop()
