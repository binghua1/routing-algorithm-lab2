CXX := g++
CXXFLAGS := -std=c++14 -O3 -I. -MMD -MP -g
LDFLAGS := 
SRCS := main.cpp LayerAssignment.cpp Router2D.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)
TARGET := router

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
