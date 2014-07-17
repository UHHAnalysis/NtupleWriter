#! /bin/sh

# This script is executed from SFrame/NtupleWriter to make the necessary object includes reachable

mkdir -p include; 
cd include;

if [ ! -h Ntuple_LinkDef.h ]; then ln -s ../interface/Ntuple_LinkDef.h .; fi;

if [ ! -h Electron.h ]; then ln -s ../Objects/Electron.h .; fi;
if [ ! -h FlavorParticle.h ]; then ln -s ../Objects/FlavorParticle.h .; fi;
if [ ! -h GenInfo.h ]; then ln -s ../Objects/GenInfo.h .; fi;
if [ ! -h GenParticle.h ]; then ln -s ../Objects/GenParticle.h .; fi;
if [ ! -h GenJetWithParts.h ]; then ln -s ../Objects/GenJetWithParts.h .; fi;
if [ ! -h GenTopJet.h ]; then ln -s ../Objects/GenTopJet.h .; fi;
if [ ! -h Jet.h ]; then ln -s ../Objects/Jet.h .; fi;
if [ ! -h MET.h ]; then ln -s ../Objects/MET.h .; fi;
if [ ! -h Muon.h ]; then ln -s ../Objects/Muon.h .; fi;
if [ ! -h Tau.h ]; then ln -s ../Objects/Tau.h .; fi;
if [ ! -h TopJet.h ]; then ln -s ../Objects/TopJet.h .; fi;
if [ ! -h Particle.h ]; then ln -s ../Objects/Particle.h .; fi;
if [ ! -h PFParticle.h ]; then ln -s ../Objects/PFParticle.h .; fi;
if [ ! -h Photon.h ]; then ln -s ../Objects/Photon.h .; fi;
if [ ! -h PrimaryVertex.h ]; then ln -s ../Objects/PrimaryVertex.h .; fi;

if [ ! -h JetProps.h ]; then ln -s ../interface/JetProps.h .; fi;
if [ ! -h Njettiness.h ]; then ln -s ../interface/Njettiness.h .; fi;
if [ ! -h Njettiness.icc ]; then ln -s ../interface/Njettiness.icc .; fi;
if [ ! -h NjettinessPlugin.h ]; then ln -s ../interface/NjettinessPlugin.h .; fi;
if [ ! -h Nsubjettiness.h ]; then ln -s ../interface/Nsubjettiness.h .; fi;
if [ ! -h Qjets.h ]; then ln -s ../interface/Qjets.h .; fi;
if [ ! -h QjetsPlugin.h ]; then ln -s ../interface/QjetsPlugin.h .; fi;

# to allow CMSSW-like includes such as "UHHAnalysis/NtupleWriter/interface/..." from within standalone sframe,
# emulate the CMSSW directory structure by creating corresponding symlinks in $SFRAME_DIR:
cd $SFRAME_DIR
mkdir UHHAnalysis
ln -s ../NtupleWriter UHHAnalysis
