#ref: http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash

#!/bin/sh

#f == all files having .png exentsion.  The following for-loop cycle through each .png file and changes it into a .jpg file 
for f in *.png; do 
    f2=$(basename $f)
    filename=${f2%.*}
    convert -quality 100 $filename.png $filename.jpg;
done