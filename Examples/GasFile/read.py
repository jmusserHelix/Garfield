import ROOT
import Garfield
import os

path = os.getenv('GARFIELD_INSTALL')

# Set up the gas.
gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('ar_80_co2_20_2T.gas')
gas.LoadIonMobility(path + '/share/Garfield/Data/IonMobility_Ar+_Ar.txt')
gas.PrintGas()

view = ROOT.Garfield.ViewMedium()
view.SetMedium(gas)
view.SetMagneticField(2.)
  
cV = ROOT.TCanvas('cV', '', 600, 600)
view.SetCanvas(cV)
view.PlotElectronVelocity()

cD = ROOT.TCanvas('cD', '', 600, 600)
view.SetCanvas(cD)
view.PlotElectronDiffusion()

cT = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cT)
view.PlotElectronTownsend()

cA = ROOT.TCanvas('cA', '', 600, 600)
view.SetCanvas(cA)
view.PlotElectronAttachment()

