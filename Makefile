CXX = g++
CXXFLAGS = -Wall -Iinclude -O3

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/file_parser.cpp $(SRC_DIR)/quad.cpp
OBJS = $(BUILD_DIR)/main.o $(BUILD_DIR)/file_parser.o $(BUILD_DIR)/quad.o
TARGET = fem.exe

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(BUILD_DIR) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)