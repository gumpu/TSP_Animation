# Usage
#
#    make .anim
#    make

sa.mp4 : 
	ffmpeg -r 10 -i %5d.png -vb 20M sa.mp4

.anim : anim.py
	python anim.py


clean:
	-rm -f *.png
	-rm -f *.mp4
