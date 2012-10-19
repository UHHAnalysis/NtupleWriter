#!/usr/bin/env python

#usage: ./myLumiCalc.py JSONname=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-191859_8TeV_PromptReco_Collisions12_JSON.txt HLTpath="HLT_PFJet320_*" outfilename=Out.root


import os,sys,time
import coral

from ROOT import gROOT, TFile, TTree, AddressOf, TString, std
from ctypes import *

gROOT.Reset()


gROOT.ProcessLine(\
       "struct MyStruct{\
       UInt_t runnr;\
       UInt_t luminosityBlock;\
       Double_t intgRecLumi;\
       Double_t HLTpresc;\
       Double_t L1presc;\
       };")
from ROOT import MyStruct



from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('JSONname',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'name of JSON file')

options.register ('HLTpath',
                  'HLT_*',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'name of HLT trigger path')

options.register ('outfilename',
                  'OutFile.root',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'name of output root file')



options.parseArguments()

os.system("lumiCalc2.py -i "+options.JSONname+" lumibyls -b stable --hltpath "+options.HLTpath+" -o TempOut.csv")
    
ofile = TFile(options.outfilename,"RECREATE")
tr = TTree("AnalysisTree","AnalysisTree")
    
s = MyStruct()

hltpath = TString("")
    

tr.Branch('run',AddressOf(s,'runnr'),'runnr/i')
tr.Branch('luminosityBlock',AddressOf(s,'luminosityBlock'),'luminosityBlock/i')
tr.Branch("HLTpath","TString",hltpath)
#tr.Branch('L1bit','TString',l1bit)
tr.Branch('intgRecLumi',AddressOf(s,'intgRecLumi'),'intgRecLumi/D')   
tr.Branch('HLTpresc',AddressOf(s,'HLTpresc'),'HLTpresc/D')
tr.Branch('L1presc',AddressOf(s,'L1presc'),'L1presc/D')

#'''
#input:  {run:[lumilsnum(0),cmslsnum(1),timestamp(2),beamstatus(3),beamenergy(4),deliveredlumi(5),recordedlumi(6),calibratedlumierror(7),{hltpath:[l1name,l1prescale,hltprescale,efflumi]},bxdata,beamdata]}
#'''
#result=[]#[run,ls,hltpath,l1bitname,hltpresc,l1presc,efflumi]


infile = open("TempOut.csv", "r")

infile.readline()
for line in infile:
    #print line
    run = line[0:line.find(":")]
    s.runnr = int(run)
    line = line[len(run)+1:]
    fill = line[0:line.find(",")]
    line = line[len(fill)+1:]
    lb = line[0:line.find(":")]
    s.luminosityBlock = int(lb)
    line = line[2*len(lb)+2:]
    hltpathname = line[0:line.find(",")]
    hltpath = TString(hltpathname)
    hltpath = TString(hltpathname)
    line = line[len(hltpathname)+1:]
    if "HLT" not in hltpathname:
        print "strange trigger name ("+hltpathname+") found for run "+run+" and lumiblock "+lb
        continue
    L1pathname =  line[0:line.find(",")]
    line = line[len(L1pathname)+1:]
    hltprescale = line[0:line.find(",")]
    s.HLTpresc = int(hltprescale)
    line = line[len(hltprescale)+1:]
    L1prescale = line[0:line.find(",")]
    s.L1presc = int(L1prescale)
    line = line[len(L1prescale)+1:]
    intLumi = line[0:line.find(",")]
    s.intgRecLumi = float(intLumi)
    #print run +"  "+lb+"  "+hltpathname 
    tr.Fill()
    
ofile.Write()
ofile.Close()  


os.system("rm TempOut.csv")

    


