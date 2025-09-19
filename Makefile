
# Compiler
CXX = clang++
CXXFLAGS = -std=c++20 -Wall -Wextra -I. -I./dependencies/eigen-3.4.0

# Directories
BUILD_DIR = build
BIN_DIR   = bin

# Source files
SRCS = main.cpp
OBJS = $(BUILD_DIR)/main.o

# Output
TARGET = $(BIN_DIR)/euclid

# === Rules ===

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

$(BUILD_DIR)/main.o: main.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean