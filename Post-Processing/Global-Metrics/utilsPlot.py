import matplotlib.pyplot as plt
import numpy as np

def PlotMetricsAgainstNVessels(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(NVessels, NTerms, NVessels, Dlim/Dlim.max(),
             NVessels, Cost/Cost.max(), NVessels, ICD/ICD.max(),
             NVessels, VAD/VAD.max(), NVessels, VSD/VSD.max(),
             NVessels, VPI/VPI.max(), NVessels, VCI/VCI.max(),
             NVessels, VDI/VDI.max())
    
    plt.legend(['NTerms', 'DLim', 'Cost', 'ICD', 'VAD', 'VSD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('Number of vessels')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
    
    plt.show()

    
def PlotMetricsAgainstDLim(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)

    plt.plot(Dlim, NTerms/NTerms.max(),
             Dlim, NVessels/NVessels.max(),
             Dlim, Cost/Cost.max(),
             Dlim, ICD/ICD.max(),
             Dlim, VAD/VAD.max(),
             Dlim, VSD/VSD.max(),
             Dlim, VPI/VPI.max(),
             Dlim, VCI/VCI.max(),
             Dlim, VDI/VDI.max())
    
    plt.legend(['NTerms', 'NVessels', 'Cost', 'ICD', 'VAD', 'VSD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('DLim')
    
    plt.gca().invert_xaxis()
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
    
    plt.show()
    
def PlotMetricsAgainstNTerms(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(NTerms, Dlim/Dlim.max(),
             NTerms, NVessels/NVessels.max(),
             NTerms, Cost/Cost.max(),
             NTerms, ICD/ICD.max(),
             NTerms, VAD/VAD.max(),
             NTerms, VSD/VSD.max(),
             NTerms, VPI/VPI.max(),
             NTerms, VCI/VCI.max(),
             NTerms, VDI/VDI.max())
    
    plt.legend(['Dlim', 'NVessels', 'Cost', 'ICD', 'VAD', 'VSD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('Number of terminal vessels')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
    
    plt.show()

def PlotMetricsAgainstGrowthStage(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(range(NVessels.size), NVessels/NVessels.max(),
             range(NVessels.size), NTerms/NTerms.max(),
             range(NVessels.size), Dlim/Dlim.max(),
             range(NVessels.size), Cost/Cost.max(),
             range(NVessels.size), ICD/ICD.max(),
             range(NVessels.size), VAD/VAD.max(),
             range(NVessels.size), VSD/VSD.max(),
             range(NVessels.size), VPI/VPI.max(),
             range(NVessels.size), VCI/VCI.max(),
             range(NVessels.size), VDI/VDI.max())

    plt.legend(['NVessels', 'DLim', 'Cost', 'ICD', 'VAD', 'VSD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('Growth stage')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
    
    plt.show()

def PlotMetricsAgainstICD(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)

    plt.plot(ICD, NVessels/NVessels.max(),
             ICD, NTerms/NTerms.max(),
             ICD, Dlim/Dlim.max(),
             ICD, Cost/Cost.max(),
             ICD, VAD/VAD.max(),
             ICD, VSD/VSD.max(),
             ICD, VPI/VPI.max(),
             ICD, VCI/VCI.max(),
             ICD, VDI/VDI.max())
    
    plt.legend(['NVessels', 'NTerms', 'DLim', 'Cost', 'VAD', 'VSD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('ICD')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
    
    plt.show()
    
def PlotMetricsAgainstVAD(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(VAD, NVessels/NVessels.max(),
             VAD, NTerms/NTerms.max(),
             VAD, Dlim/Dlim.max(),
             VAD, Cost/Cost.max(),
             VAD, ICD/ICD.max(),
             VAD, VSD/VSD.max(),
             VAD, VPI/VPI.max(),
             VAD, VCI/VCI.max(),
             VAD, VDI/VDI.max())
    
    plt.legend(['NVessels', 'NTerms', 'DLim', 'Cost', 'ICD', 'VSD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('VAD')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
        
    plt.show()
    
def PlotMetricsAgainstVSD(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(VSD, NVessels/NVessels.max(),
             VSD, NTerms/NTerms.max(),
             VSD, Dlim/Dlim.max(),
             VSD, Cost/Cost.max(),
             VSD, ICD/ICD.max(),
             VSD, VAD/VAD.max(),
             VSD, VPI/VPI.max(),
             VSD, VCI/VCI.max(),
             VSD, VDI/VDI.max())
    
    plt.legend(['NVessels', 'NTerms', 'DLim', 'Cost', 'ICD', 'VAD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('VSD')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
        
    plt.show()
    
def PlotMetricsAgainstVPI(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(VPI, NVessels/NVessels.max(),
             VPI, NTerms/NTerms.max(),
             VPI, Dlim/Dlim.max(),
             VPI, Cost/Cost.max(),
             VPI, ICD/ICD.max(),
             VPI, VAD/VAD.max(),
             VPI, VSD/VSD.max(),
             VPI, VCI/VCI.max(),
             VPI, VDI/VDI.max())
    
    plt.legend(['NVessels', 'NTerms', 'DLim', 'Cost', 'ICD', 'VAD', 'VPI', 'VCI', 'VDI'])
    plt.xlabel('VPI')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
        
    plt.show()
    
def PlotMetricsAgainstVCI(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(VCI, NVessels/NVessels.max(),
             VCI, NTerms/NTerms.max(),
             VCI, Dlim/Dlim.max(),
             VCI, Cost/Cost.max(),
             VCI, ICD/ICD.max(),
             VCI, VAD/VAD.max(),
             VCI, VSD/VSD.max(),
             VCI, VPI/VPI.max(),
             VCI, VDI/VDI.max())
    
    plt.legend(['NVessels', 'NTerms', 'DLim', 'Cost', 'ICD', 'VAD', 'VSD', 'VCI', 'VDI'])
    plt.xlabel('VCI')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
    
    plt.show()

def PlotMetricsAgainstVDI(metricsFile, logscale = True):
    NTerms, Dlim, Cost, NVessels, ICD, VAD, VSD, VPI, VCI, VDI = np.loadtxt(metricsFile, unpack=True)
    
    plt.plot(VDI, NVessels/NVessels.max(),
             VDI, NTerms/NTerms.max(),
             VDI, Dlim/Dlim.max(),
             VDI, Cost/Cost.max(),
             VDI, ICD/ICD.max(),
             VDI, VAD/VAD.max(),
             VDI, VSD/VSD.max(),
             VDI, VPI/VPI.max(),
             VDI, VCI/VCI.max())
    
    plt.legend(['NVessels', 'NTerms', 'DLim', 'Cost', 'ICD', 'VAD', 'VSD', 'VPI', 'VCI'])
    plt.xlabel('VDI')
    plt.ylabel('Normalized metrics')
    if logscale:
        plt.yscale('log')
        
    plt.show()
