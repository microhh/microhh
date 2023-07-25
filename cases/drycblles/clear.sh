rm u.*
rm v.*
rm w.*
rm th*.*
rm time.*
rm d*d*.*
rm grid.*
rm fftwplan.*
#rm drycblles_input.nc
find . -name "drycblles.*" -not -name "*.ini" -exec rm {} \;
