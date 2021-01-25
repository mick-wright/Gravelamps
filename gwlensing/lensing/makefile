CC = g++
TARGET = plCalc
OBJECTS = pointlens.cxx
CPPFLAGS = -O3 -larb -lflint

$(TARGET): $(OBJECTS)
	@echo "Building point lensing code"
	$(CC) $(OBJECTS) -o $(TARGET) $(CPPFLAGS)
	@cp $(TARGET) /usr/local/bin/$(TARGET) 

clean:
	@rm $(TARGET)
	@rm /usr/local/bin/$(TARGET) 

