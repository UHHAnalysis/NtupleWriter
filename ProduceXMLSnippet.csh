#!/bin/csh

set outfile = `basename $1`'.xml'
echo "Creating XML snippet: "$outfile
rm $outfile
foreach file (`ls /pnfs/cms/WAX/11/store/user/b2g12006/$1`)
    echo '<In FileName="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/b2g12006/'$1'/'$file'" Lumi="0.0"/>' >> $outfile
end
