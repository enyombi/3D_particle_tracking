#ref: http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash

#!/bin/sh

#f == all files having .jpg exentsion.  The following for-loop cycle through each .jpg file, adds 1 to the interative variable 'd' and changes the filenamae to a 3 digit number having the value of 'd' (e.g., when d==1, mv import_p01.jpg 001.jpg; when d==2, mv import_p02.jpg 002.jpg, etc.) 

d=0
for f in *.jpg; do 
    f2=$(basename $f)
    filename=${f2%.*}
    d=$((d+1))
    filenameNew=$(printf "%03d" $d)
    mv "$filename.jpg" "$filenameNew.jpg"
done