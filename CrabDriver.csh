#!/bin/csh
if ($#argv == 0) then
    echo "Usage: source CrabDriver.csh {DataSetPath} {DataSetName} [LumiMask]"
    echo " {..} are required, LumiMask] is required for data samples"
    echo " Example: source CrabDiver.csh /ElectronHad/pturner-Run2012A-13Jul2012-v1_TLBSM_53x_v2-e3fb55b810dc7a0811f4c66dfa2267c9/USER ElectronHad_Run2012A MyJSON.txt"
    exit
endif

set file="crab_NtupleWriter_Template.cfg"
set DSPath = $1
set DataSetName = $2
set isData = 0
set LumiMask = ""
if ($#argv == 3) then
    set LumiMask = $3
    set isData = 1
endif

set OutDir = "/store/user/b2g12006/53xNTuples-v2/"${DataSetName}

if (`echo $DSPath | grep "StoreResults"` == `echo $DSPath`) then
    set StoreResults = 1
else
    set StoreResults = 0
endif

if ($StoreResults == 1) then
    echo "s|\#dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_prod_global/|dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_prod_global/|g" > sedscript5
else
    echo "s|\#dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/|dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/|g" > sedscript5
endif

if ($isData == 1) then
    echo "s|\#lumis_per_job|lumis_per_job|g" > sedscript6
    echo "s|\#total_number_of_lumis|total_number_of_lumis|g" >> sedscript6
    echo "s|\#lumi_mask|lumi_mask|g" >> sedscript6
else
    echo "s|\#total_number_of_events|total_number_of_events|g" > sedscript6
    echo "s|\#events_per_job|events_per_job|g" >> sedscript6
endif



echo "s|{InsertDatasetPathHere}|"${DSPath}"|g" > sedscript1
echo "s|{LumiMask}|"${LumiMask}"|g" > sedscript2
echo "s|{RemoteOutputDir}|"${OutDir}"|g" > sedscript3
echo "s|{LocalCrabDir}|"${DataSetName}"|g" > sedscript4
mkdir crabcfg
cat $file | sed -f sedscript1 > outtemp1
cat outtemp1 | sed -f sedscript2 > outtemp2
cat outtemp2 | sed -f sedscript3 > outtemp3
cat outtemp3 | sed -f sedscript4 > outtemp4
cat outtemp4 | sed -f sedscript5 > outtemp5
cat outtemp5 | sed -f sedscript6 > ${PWD}/crabcfg/crab_"${DataSetName}".cfg
rm sedscript*
rm outtemp*
