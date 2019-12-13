#!/bin/bash

S = $1
#ffmpeg -i DJI_0023.MOV -s 1920x1080 -vf "drawtext=fontfile=/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf:expansion=normal: text='frame %{n}\\: pts=%{pts \\: hms}': fontcolor=white:fontsize=48: x=7: y=7" -vcodec libx264 -vb 2000k -strict -2 -preset ultrafast -f mp4 output.mp4

ffmpeg -i $1 -r 30 -t 100 image-%d.jpeg 
ls -1v image-*.jpeg > images.txt
ffprobe -f lavfi -i "movie="${1}",fps=fps=30[out0]" -show_frames -show_entries frame=pkt_pts_time -of csv=p=0 > frames.txt
paste images.txt frames.txt > combined.txt