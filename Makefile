# Compilatore e flag
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
INCLUDES = -Iinclude

# File sorgenti
SRCS = src/main.cpp src/CSVParser.cpp src/VCFReconstructor.cpp

# Eseguibile finale
TARGET = build/test_parser

# Regola di default quando scrivi solo "make"
all: $(TARGET)

# Come costruire l'eseguibile
$(TARGET): $(SRCS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS) -o $(TARGET)
	@echo "Compilazione completata! Eseguibile salvato in $(TARGET)"

# Regola per ripulire i file compilati (comando: make clean)
clean:
	rm -rf build/*
	@echo "Pulizia completata."