#!/bin/bash
for i in $1/img/ForagingModel/*.png
do
    PNGNAME=$i
    echo $PNGNAME
    JPGNAME=$(echo $PNGNAME | sed -r 's/.png/.jpg/g')
    convert $PNGNAME -background white -alpha remove -alpha off $JPGNAME
done
ffmpeg -y -r 30 -f image2 -i $1/img/ForagingModel/ForagingModel-t%d0.jpg -vcodec libx264 -crf 25 -pix_fmt yuv420p $1/mp4/ForagingModel.mp4
