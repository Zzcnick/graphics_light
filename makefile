all: Picture.java
	javac Picture.java && java Picture test.scr && echo "Making tri.gif..." && convert *.ppm -loop 0 -delay 2 tri.gif && animate tri.gif

test: Picture.java
	javac Picture.java && java Picture test2.scr && echo "Making spin.gif..." &&  convert *.ppm -loop 0 -delay 20 spin.gif && display spin.gif

test2: Picture.java
	javac Picture.java && java Picture && display triz.ppm && display tri.ppm

run2: Picture.java
	make && echo "Running..." && java Picture && make png && rm out.ppm && echo "Saved as out.png" && display out.png

script:
	make && echo "Running script..." && java Picture test.scr && make png && rm out.ppm && echo "Saved as out.png" && display out.png

clean: 
	rm *.class *.ppm *~ *.jpg *.png

jpg: out.ppm
	java Picture; \
	convert out.ppm out.jpg

png: out.ppm
	java Picture; \
	convert out.ppm out.png
