echo -n 'Please enter directory where ntuples are stored:'
set dir = $<
foreach file (`ls $dir`)
   echo '<In FileName="'$PWD'/'$dir'/'$file'" Lumi="0.0"/>' >> `basename $dir`'.xml'
end
