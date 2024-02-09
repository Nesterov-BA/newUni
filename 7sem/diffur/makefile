# Variables
CC = g++
CFLAGS = -Wall -Wextra -O2 -g
OBJDIR = obj
SRCDIR = src
BINDIR = bin
TARGET = $(BINDIR)/a.out
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
DEPS = $(OBJS:.o=.d)

# Phony Targets
.PHONY: all clean

# Default Target
all: $(TARGET)

# Compile Object Files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Generate Dependencies
$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	$(CC) -MM $< > $@

# Link Object Files
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Include Dependencies
-include $(DEPS)

# Clean Target
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d $(TARGET)