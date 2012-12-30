FLAGS = -g
LIBS = -lglut -lGLU -lGL -lm
GLUI_LIBS = -lglui -lglut -lGLU -lGL -lm


all: 3dcube

3dcube: 3dcube.cpp
	g++ 3dcube.cpp -o 3dcube $(FLAGS) $(GLUI_LIBS)

test: all
	./3dcube

clean: all
	rm 3dcube
