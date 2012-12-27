FLAGS = -g
LIBS = -lglut -lGLU -lGL -lm
GLUI_LIBS = -lglui -lglut -lGLU -lGL -lm


all: hw2

hw2: hw2.cpp
	g++ hw2.cpp -o hw2 $(FLAGS) $(GLUI_LIBS)

test: all
	./hw2

clean: all
	rm -rf hw2