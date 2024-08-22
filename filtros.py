#filtros.py
import sys 
import os 
import numpy as np
import math 
from sympy import *
from PySide6.QtWidgets import *
from PySide6.QtGui import *
from PySide6.QtCore import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class Filtros(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PID Control Interface")
        self.setup_ui()

    def setup_ui(self):
        self.setMinimumSize(500, 650)

        self.main_frame = QFrame()
        self.main_layout = QVBoxLayout(self.main_frame)

        self.central_frame = QFrame()
        self.central_layout = QHBoxLayout(self.central_frame)

        self.left_frame = QFrame()
        self.left_layout = QVBoxLayout(self.left_frame)
        self.label1 = QLabel("Maximum variation in the passband: (Am)")
        self.label2 = QLabel("Frequency that delimits the pass band: ωp")
        self.label3 = QLabel("Minimum attenuation in the rejection band: (Amin)")
        self.label4 = QLabel("Frequency that delimits the rejection range: ωs")
        self.label5 = QLabel("Filter type:")
        self.label6 = QLabel("Topology type:")
        self.label7 = QLabel("Gain (k): ")
        self.label7.hide()
        self.label_empty = QLabel(" ")
        self.left_layout.addWidget(self.label1)
        self.left_layout.addWidget(self.label2)
        self.left_layout.addWidget(self.label3)
        self.left_layout.addWidget(self.label4)
        self.left_layout.addWidget(self.label5)
        self.left_layout.addWidget(self.label6)
        self.left_layout.addWidget(self.label7)
        self.left_layout.addWidget(self.label_empty)

        # create a vertical line as separator between left and right side
        self.vertical_line = QFrame()
        self.vertical_line.setFrameShape(QFrame.VLine)
        self.vertical_line.setFrameShadow(QFrame.Sunken)

        self.right_frame = QFrame()
        self.right_layout = QVBoxLayout(self.right_frame)
        self.am_input = QLineEdit("1")
        self.wp_input = QLineEdit("10000")
        self.amin_input = QLineEdit("25")
        self.ws_input = QLineEdit("15000")
        self.filter_type = QComboBox()
        self.filter_type.addItems(['Butterworth','Chebyshev'])
        self.topology_type = QComboBox()
        self.topology_type.addItems(['Sallen-Key', 'Multiple Feedback'])
        self.topology_type.currentIndexChanged.connect(self.gain_change)
        self.gain_input = QLineEdit("1")
        self.gain_input.hide()
        self.button_update = QPushButton("Update")
        self.button_update.clicked.connect(self.select_filter)
        self.right_layout.addWidget(self.am_input)
        self.right_layout.addWidget(self.wp_input)
        self.right_layout.addWidget(self.amin_input)
        self.right_layout.addWidget(self.ws_input)
        self.right_layout.addWidget(self.filter_type)
        self.right_layout.addWidget(self.topology_type)
        self.right_layout.addWidget(self.gain_input)
        self.right_layout.addWidget(self.button_update)

        self.central_layout.addWidget(self.left_frame)
        self.central_layout.addWidget(self.vertical_line)
        self.central_layout.addWidget(self.right_frame)

        self.bottom_frame = QFrame(self)
        self.bottom_layout = QHBoxLayout(self.bottom_frame)

        self.under_bottom_frame = QFrame()
        self.under_bottom_layout = QHBoxLayout(self.under_bottom_frame)
        self.button_plan_s = QPushButton("Plan s")
        self.button_plan_s.clicked.connect(self.show_plan_s)
        self.button_bode_diagram = QPushButton("Bode Diagram")
        self.button_bode_diagram.clicked.connect(self.select_bode_diagram)
        self.button_components = QPushButton("Get components")
        self.button_components.clicked.connect(self.select_topology)
        self.under_bottom_layout.addWidget(self.button_plan_s)
        self.under_bottom_layout.addWidget(self.button_bode_diagram)
        self.under_bottom_layout.addWidget(self.button_components)

        self.main_layout.addWidget(self.central_frame)
        self.main_layout.addWidget(self.bottom_frame)
        self.main_layout.addWidget(self.under_bottom_frame)
        self.setCentralWidget(self.main_frame)

    def gain_change (self):
        selected_unit = self.topology_type.currentText()
        if selected_unit == 'Multiple Feedback':
            self.label7.show()
            self.gain_input.show()
        else:
            self.label7.hide()
            self.gain_input.hide()

    def select_filter (self):
        am = self.am_input.text()
        wp = self.wp_input.text()
        amin = self.amin_input.text()
        ws = self.ws_input.text()
        if self.gain_input.isEnabled():
            gain_text = self.gain_input.text()
        else:
            gain_text = "1"
        filter = self.filter_type.currentText()
        if filter == 'Butterworth':
            self.start_butterworth(am,wp,amin,ws,gain_text)
        if filter == 'Chebyshev':
            self.start_chebyshev(am,wp,amin,ws,gain_text)

    def start_butterworth(self, am, wp, amin, ws, k):
        k = float(k)
        am = float(am)
        wp = 2*math.pi*float(wp)
        amin = float(amin)
        ws = 2*math.pi*float(ws)
        self.wp = wp
        self.ws = ws
        #print(am,wp,amin,ws)
        ep = sqrt(10**(am/10)-1)
        self.ep = float(ep)
        #print("ep = ", ep)
        n = math.log10((10**(amin/10)-1)/((ep)**2))/(2*math.log10(ws/wp))
        self.n = n
        #print("n = ",n)
        nint = 0
        while nint <n:
            nint+=1
        #print("nint = ",nint)
        w = wp*2*math.pi*((1/ep)**(1/nint))
        self.polos = np.zeros(nint, dtype=complex)
        j = len(self.polos)-1
        if nint % 2 == 0:
            theta = math.pi/(2*nint)
            #print("theta = ", theta*180/math.pi)
            theta2 = math.pi/nint
            #print("thetha2 = ", theta2*180/math.pi)
            for i in range(0, nint//2):
                theta3 = math.pi/2-(theta+theta2*i)
                self.polos[i] = complex(w*-cos(theta3),w*sin(theta3))
                self.polos[j-i] = complex(w*-cos(theta3),w*-sin(theta3))
        else:
            theta = math.pi/(2*nint)
            #print("theta = ", theta*180/math.pi)
            theta2 = math.pi/nint
            #print("thetha2 = ", theta2*180/math.pi)
            for i in range(0, nint//2):
                #print('i = ',i,'; j-i = ',j-i)
                theta3 = math.pi/2-(theta+theta2*i)
                self.polos[i] = complex(w*-cos(theta3),w*sin(theta3))
                self.polos[j-i] = complex(w*-cos(theta3),w*-sin(theta3))
            self.polos[nint//2] = complex(-w,0)
        self.show_ft_butteworth(self.polos, k)
        #print("polos = ",self.polos)

    def show_ft_butteworth(self,polos,k):
        self.polos = polos
        j = len(self.polos)
        if j%2==0:
            kpart = k**(1/(j//2))
        else:
            kpart = k**(1/((j+1)//2))
        equation_parts = []
        for i in range(0, j//2):
            p1 = polos[i]
            w = abs(p1)
            zeta = -p1.real/w
            num = kpart * w**2
            denom = f"s^2 + {2*zeta*w:.2f} s + {w**2:.2f}"
            equation_parts.append(rf"\frac{{{num:.2f}}}{{{denom}}}")
        if j%2==1:
            p_extra = polos[j//2]
            num = kpart
            denom = f"s + {-p_extra.real:.2f}"
            equation_parts.append(rf"\frac{{{num}}}{{{denom:}}}")

        latex = " * ".join(equation_parts)
        latex = f"H(s) = {latex}"
        latex = f"${latex}$"
        #print(latex)
        fig = Figure(figsize=(9, 3))
        ax = fig.add_subplot(111, facecolor='#404040')
        ax.text(0.5, 0.5,latex, fontsize=15, ha='center', va='center', color='white')
        ax.axis('off')
        fig.patch.set_facecolor('#404040')
        canvas = FigureCanvas(fig)

        for i in reversed(range(self.bottom_layout.count())):
            widget_to_remove = self.bottom_layout.itemAt(i).widget()
            self.bottom_layout.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)

        self.bottom_layout.addWidget(canvas)
        self.bottom_frame.updateGeometry()
        self.bottom_frame.repaint()

    def start_chebyshev(self, am, wp, amin, ws, gain):
        am = float(am)
        wp = 2*math.pi*float(wp)
        amin = float(amin)
        ws = 2*math.pi*float(ws)
        gain = float(gain)
        self.wp = wp
        self.ws = ws
        ep = sqrt(10**(am/10)-1)
        n = math.acosh(sqrt((10**(amin/10)-1)/ep**2)/math.acosh(ws/wp))
        self.n = n
        self.ep = float(ep)
        #print("ep = ", ep)
        #print("n = ",n)
        nint = 0
        while nint <n:
            nint+=1
        #print("nint = ",nint)
        self.polos = np.zeros(nint, dtype=complex)
        for k in range(1,nint+1):
           a=math.sin(((2*k-1)/n)*math.pi/2)
           b=math.sinh((1/n)*math.asinh(1/ep))
           c=math.cos(((2*k-1)/n)*math.pi/2)
           d=math.cosh((1/n)*math.asinh(1/ep))
           self.polos[k-1] = complex(-wp*a*b,wp*c*d)
        self.show_ft_chebyshev(self.polos,n,ep,wp,gain)
        
    def show_ft_chebyshev(self,polos,n,ep,wp,k):
        j = len(polos)
        equation_parts = []
        num = k * wp**n
        denom_parts = [f"{ep *(2**(n-1)):.2f}"]
        for i in range(0, j):
            p1 = -polos[i]
            denom_parts.append(f" * (s + {p1:.2f}) ")
        denom_str = " ".join(denom_parts)

        equation_parts.append(rf"\frac{{{num:.2f}}}{{{denom_str}}}")
        latex = " * ".join(equation_parts)
        latex = f"H(s) = {latex}"
        latex = f"${latex}$"
        #print(latex)
        fig = Figure(figsize=(9, 3))
        ax = fig.add_subplot(111, facecolor='#404040')
        ax.text(0.5, 0.5,latex, fontsize=15, ha='center', va='center', color='white')
        ax.axis('off')
        fig.patch.set_facecolor('#404040')
        canvas = FigureCanvas(fig)

        for i in reversed(range(self.bottom_layout.count())):
            widget_to_remove = self.bottom_layout.itemAt(i).widget()
            self.bottom_layout.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)

        self.bottom_layout.addWidget(canvas)
        self.bottom_frame.updateGeometry()
        self.bottom_frame.repaint()

    def show_plan_s(self):

        # Configurar o gráfico
        plt.figure(figsize=(8, 8))
        plt.axhline(0, color='black', lw=0.5)  # Eixo x
        plt.axvline(0, color='black', lw=0.5)  # Eixo y

        plt.scatter(self.polos.real, self.polos.imag, marker='x', color='red', label='Polos')

        # Adicionar etiquetas e grade
        plt.xlabel('Parte Real')
        plt.ylabel('Parte Imaginária')
        plt.title('Diagrama de Polos')
        plt.legend()
        plt.grid(True)

        plt.show()

    def select_bode_diagram(self):

        filter = self.filter_type.currentText()
        if filter == 'Butterworth':
            self.show_bode_diagram_butterworth()
        if filter == 'Chebyshev':
            self.show_bode_diagram_chebyshev()

    def show_bode_diagram_butterworth(self):
        
        if not hasattr(self, 'ep') or not hasattr(self, 'wp') or not hasattr(self, 'n'):
            print("Os parâmetros necessários não estão definidos.")
            return

        w = np.logspace(2, 10, 1000)
        try:
            tjw = 1 / np.sqrt(1 +(self.ep**2)*((w/self.wp)**(2*self.n)))
            tjw_db = 20 * np.log10(tjw)
        except Exception as e:
            print(f"Erro no cálculo: {e}")
            return

        #construindo o segundo diagrama de fase
        try:
            j = len(self.polos)
            fiw = np.zeros_like(w)
            for polo in self.polos:
                wp = abs(polo)
                fiw -= np.arctan(w/wp)
            fiw_deg = np.degrees(fiw)
            fiw_deg = (fiw_deg + 180) % 360 - 180
        except Exception as e:
            print(f"Erro no cálculo da fase: {e}")
            return

        # Plotar o diagrama de Bode (magnitude)
        plt.figure(figsize=(18, 8))

        plt.subplot(1,2,1)
        plt.semilogx(w, tjw_db)
        plt.xlabel('Frequência (rad/s)')
        plt.ylabel('Magnitude (dB)')
        plt.title('Diagrama de Bode - Magnitude')
        plt.grid(which='both', linestyle='--', linewidth=0.5)
        
        plt.subplot(1,2,2)
        plt.semilogx(w, fiw_deg)
        plt.xlabel('Frequência (rad/s)')
        plt.ylabel('Fase (Graus)')
        plt.title('Diagrama de Fase')
        plt.grid(which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout()
        plt.show()



    def show_bode_diagram_chebyshev(self):
        if not hasattr(self, 'ep') or not hasattr(self, 'wp') or not hasattr(self, 'n'):
            print("Os parâmetros necessários não estão definidos.")
            return

        w = np.logspace(-1, 6, 1000)
        try:
            tjw = np.zeros_like(w)
            for i in range(len(w)):
                if (w[i]>=self.wp):
                    tjw[i] = 1 / np.sqrt(1+(self.ep**2)*np.cosh((self.n*np.arccosh(w[i]/self.wp)))**2)
                else:
                    tjw[i] = 1 / np.sqrt(1+(self.ep**2)*np.cos((self.n*np.arccos(w[i]/self.wp)))**2)
            tjw_db = 20 * np.log10(tjw)
        except Exception as e:
            print(f"Erro no cálculo: {e}")
            return

        #diagrama de fase
        try:
            j = len(self.polos)
            fiw = np.zeros_like(w)
            for polo in self.polos:
                wp = abs(polo)
                fiw -= np.arctan(w/wp)
            fiw_deg = np.degrees(fiw)
            fiw_deg = (fiw_deg + 180) % 360 - 180
        except Exception as e:
            print(f"Erro no cálculo da fase: {e}")
            return


        # Plotar o diagrama de Bode (magnitude)
        plt.figure(figsize=(18, 8))

        plt.subplot(1,2,1)
        plt.semilogx(w, tjw_db)
        plt.xlabel('Frequência (rad/s)')
        plt.ylabel('Magnitude (dB)')
        plt.title('Diagrama de Bode - Magnitude')
        plt.grid(which='both', linestyle='--', linewidth=0.5)
        plt.show()

        plt.subplot(1,2,2)
        plt.semilogx(w, fiw_deg)
        plt.xlabel('Frequência (rad/s)')
        plt.ylabel('Fase (Graus)')
        plt.title('Diagrama de Fase')
        plt.grid(which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout()
        plt.show()

    def select_topology (self):
        topology = self.topology_type.currentText()
        gain = self.gain_input.text()
        if topology == 'Sallen-Key':
            self.get_components_sallenkey()
        if topology == 'Multiple Feedback':
            self.get_components_multiple_feedback(gain)

    def get_components_sallenkey (self):
        k = len(self.polos)
        r1_parts = []
        r2_parts = []
        c1_parts = []
        c2_parts = []
        for i in range (0,k//2):
            p1 = self.polos[i].real
            p2 = self.polos[i].imag
            ksi= abs(p2)/abs(p1)
            w0 = abs(self.polos[i])
            Q = 1/(2*ksi)
            n = 4*Q**2 + 1
            m1 = (-(2-n/Q**2)+sqrt((2-n/Q**2)**2-4))/2
            m2 = (-(2-n/Q**2)-sqrt((2-n/Q**2)**2-4))/2
            c1 = 10*10**-9
            c2 = n*c1
            r2 = 1/(c1*w0*sqrt(m1*n))
            r1 = m1*r2
            r1_parts.append(f"{r1:.2f} ")
            r2_parts.append(f"{r2:.2f} ")
            c1_parts.append(f"{c1:.2e} ")
            c2_parts.append(f"{c2:.2e} ") 
        if k%2==1:
            p1 = -self.polos[k//2].real
            c1 = 10*10**-9
            r1 = 1 / (c1*p1)
            r1_parts.append(f"{r1:.2f} ")
            c1_parts.append(f"{c1:.2e} ")
        r1_text = f'; '.join(r1_parts)
        r2_text = f'; '.join(r2_parts)
        c1_text = f'; '.join(c1_parts)
        c2_text = f'; '.join(c2_parts)
        '''
        print (r1_text)
        print (r2_text)
        print (c1_text)
        print (c2_text)'''
        topology = 'sallenkey'
        plot_window = PlotWindow()
        plot_window.writelabels(r1_text,r2_text,c1_text,c2_text,None,topology)
        plot_window.exec()

    def get_components_multiple_feedback(self,k):
        l = len(self.polos)
        k = float(k)
        r1_parts = []
        r2_parts = []
        r3_parts = []
        c1_parts = []
        c2_parts = []
        kpart = k**(1/(l//2))
        for i in range (0,l//2):
            p1 = self.polos[i].real
            p2 = self.polos[i].imag
            ksi = abs(p2)/abs(p1)
            w0 = abs(self.polos[i])
            c1 = 10*10**-9
            g = 1+kpart
            q = 1/(2*ksi)
            n = 4*g*q**2
            a= (2/g)-(n/((g**2)*(q**2)))
            m = -a+sqrt((a**2)-4/g**2)
            r2 = 1/(w0*c1*sqrt(m*n))
            r3 = m*r2
            r1 = r2/kpart 
            c2 = n*c1 
            r1_parts.append(f"{r1:.2f}; ")
            r2_parts.append(f"{r2:.2f}; ")
            r3_parts.append(f"{r3:.2f}; ")
            c1_parts.append(f"{c1:.2e}; ")
            c2_parts.append(f"{c2:.2e}; ")
        if k%2==1:
            p1 = -self.polos[k//2].real
            c1 = 10*10**-9
            r1 = 1 / (c1*p1)
            r1_parts.append(f"{r1:.2f} ")
            c1_parts.append(f"{c1:.2e} ")

        r1_text = f'; '.join(r1_parts)
        r2_text = f'; '.join(r2_parts)
        r3_text = f'; '.join(r3_parts)
        c1_text = f'; '.join(c1_parts)
        c2_text = f'; '.join(c2_parts)

        topology = 'multiple_feedback'
        plot_window = PlotWindow()
        plot_window.writelabels(r1_text,r2_text,c1_text,c2_text,r3_text,topology)
        plot_window.exec()

class PlotWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("RC parameters")
        self.setup_ui()

    def setup_ui(self):

        self.main_frame = QFrame()
        self.main_layout = QVBoxLayout(self.main_frame)

        self.label_frame = QFrame()
        self.label_layout = QVBoxLayout(self.label_frame)
        self.label_r1 = QLabel("R1: ")
        self.label_r2 = QLabel("R2: ")
        self.label_r3 = QLabel("R3: ")
        self.label_r3.hide()
        self.label_c1 = QLabel("C1: ")
        self.label_c2 = QLabel("C2: ")
        self.label_layout.addWidget(self.label_r1)
        self.label_layout.addWidget(self.label_r2)
        self.label_layout.addWidget(self.label_r3)
        self.label_layout.addWidget(self.label_c1)
        self.label_layout.addWidget(self.label_c2)

        self.bottom_frame = QFrame()
        self.bottom_layout = QHBoxLayout(self.bottom_frame)
        self.back_button = QPushButton("Back")
        self.back_button.clicked.connect(self.go_back)
        self.bottom_layout.addWidget(self.back_button)

        self.main_layout.addWidget(self.label_frame)
        self.main_layout.addWidget(self.bottom_frame)

        self.setLayout(self.main_layout)

    def writelabels(self,r1,r2,c1,c2,r3,topology):

        self.label_r1.setText(f"R1: {r1} Ω ")
        self.label_r2.setText(f"R2: {r2} Ω ")
        self.label_c1.setText(f"C1: {c1} F ")
        self.label_c2.setText(f"C2: {c2} F ") 

        if topology == 'multiple_feedback':
            self.label_r3.setText(f"R3: {r3} Ω")
            self.label_r3.show()
        else:
            self.label_r3.hide()
            
    def go_back(self):
        self.close()

def apply_stylesheet(app, stylesheet_path):
    print(f"Attempting to apply stylesheet from: {stylesheet_path}")
    try:
        with open(stylesheet_path, "r") as f:
            stylesheet = f.read()
        app.setStyleSheet(stylesheet)
        print("Stylesheet applied successfully.")
    except FileNotFoundError:
        print(f"Error: The file '{stylesheet_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def create_and_show_window():
    app = QApplication(sys.argv)  # Create Aplicacion
    apply_stylesheet(app,"style.qss") #Apply the css
    window = Filtros()  # Create Windows
    window.show()  # Show Windows
    sys.exit(app.exec())  # Execute application loop 

def start_application():
    create_and_show_window()

if __name__ == "__main__":
    start_application()
