#######################################################################################################

# Mac OS X
#INCLUDE_PATH      = -I/usr/local/include/ -I/usr/local/include/eigen3/
#LIBRARY_PATH      = -L/usr/local/lib/
#OPENGL_LIBS       = -framework OpenGL -framework GLUT

# # Linux
INCLUDE_PATH      = -I/usr/local/include/ -I/usr/local/include/eigen3/
LIBRARY_PATH      = -L/usr/local/lib/
OPENGL_LIBS       = -lglut -lGL -lX11

# # Windows / Cygwin
# INCLUDE_PATH      = -I/usr/include/opengl
# LIBRARY_PATH      = -L/usr/lib/w32api
# OPENGL_LIBS       = -lglut32 -lopengl32

#######################################################################################################

TARGET = stf
CC = g++
LD = g++
CFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated -pedantic $(INCLUDE_PATH) -I./include -I./src -DNDEBUG
LFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated -Werror -pedantic $(LIBRARY_PATH) -DNDEBUG
LIBS = $(OPENGL_LIBS)

OBJS = main.o
HEADERS = BodyForce.h Boundary.h Environment.h MaterialModel.h Node.h Particle.h Element.h MaterialPoint.h Objects.h Solid.h Damping.h STF.h
WEBOBJS = index.js index.wasm index.html

default: $(TARGET)

.PHONY: web

all: clean $(TARGET)

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LIBS) -o $(TARGET)

web:
	emcc main.cpp -s WASM=1 -s LEGACY_GL_EMULATION=1 -lglut -lGLU -lGL -o index.html

main.o: main.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c main.cpp -o main.o

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
	rm -f $(TARGET).exe
	rm -f $(WEBOBJS)

