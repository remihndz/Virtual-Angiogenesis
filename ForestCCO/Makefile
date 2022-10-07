# Variable Definitions
CC = g++
FLAGS = -lm --std=c++17
SRC = ~/ForestCCO/src
SRC_REMI = ./Geometry/src
INCLUDE_REMI = ./Geometry/include
EXTENSION = cc
EXEC_MACULA = macula
BASE_FILES = $(SRC)/progress/*.$(EXTENSION) $(SRC)/geometry/*.$(EXTENSION) $(SRC)/domain/*.$(EXTENSION) $(SRC)/tree/*.$(EXTENSION) $(SRC)/cco/*.$(EXTENSION) $(SRC)/morphometry/*.$(EXTENSION) $(SRC_REMI)/*.$(EXTENSION)
FILES_MACULA = $(EXEC_MACULA).$(EXTENSION) $(BASE_FILES) $(SRC)/forest/*.$(EXTENSION)
INCLUDES = -I $(SRC) -I $(SRC)/progress -I $(SRC)/forest -I $(SRC)/cco -I $(SRC)/cco/interface -I $(SRC)/tree -I $(SRC)/domain -I $(SRC)/morphometry -I $(SRC)/voronoi -I $(INCLUDE_REMI)

# Compiling rules.
all: macula

# Compiling Macula 
macula:
	@echo "Compiling $(EXEC_MACULA)..."
	$(CC) -o $(EXEC_MACULA) $(FILES_MACULA) $(INCLUDES) $(FLAGS)
	@echo ""

# Clean binaries
clean:
	@rm -f $(EXEC_MACULA)
	@echo "All binaries cleaned up!"
