import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def read_data(f):
    data = []
    file = open(f, 'r')
    for line in file:
        line_list=line.split(",")
        data.extend(line_list)
    data.pop()   
    file.close()
    return data
def create_list(data):
    x =[]
    y =[]
    for dp in data:
        x_y = []
        x_y = dp.split(":")
        x.append(float(x_y[0]))
        y.append(float(x_y[1]))
    return x , y
        
def plot_data(x1,y1,x2,y2):
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Harmonic Oscillator Basis')
    plt.plot(x1,y1, linewidth=1.0)
    plt.plot(x2,y2, linewidth=1.0)
    plt.axis([x1.min(), x1.max(), y1.min(), y1.max()])
    
    pp = PdfPages("plot.pdf")
    plt.savefig(pp, format="pdf")
    pp.close()
    
    plt.show()
    
def main():
    data = read_data("bis1.txt")
    x_,y_ = create_list(data)
    x1 = np.asarray(x_)
    y1 = np.asarray(y_)
    
    data = read_data("bis2.txt")
    x_,y_ = create_list(data)
    x2 = np.asarray(x_)
    y2 = np.asarray(y_)
    
    plot_data(x1,y1,x2,y2)
    
main()                                                                                                                                                                                                                                 