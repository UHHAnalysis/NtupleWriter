#!/bin/csh

echo -n 'Please enter directory where ntuples are stored:/pnfs/cms/WAX/11/store/user/b2g12006/53xNTuples-v2/'
set dir = $<
set outfile = `basename $dir`'.xml'
echo "Creating XML snippet: "$outfile
rm $outfile
foreach file (`ls /pnfs/cms/WAX/11/store/user/b2g12006/$dir`)
    echo '<In FileName="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/b2g12006/'$dir'/'$file'" Lumi="0.0"/>' >> $outfile
end
