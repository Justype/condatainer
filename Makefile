.PHONY: build clean test

# Binary name
BINARY_NAME=condatainer_go
BINARY_PATH=bin/$(BINARY_NAME)

# Build the binary
build:
	@echo "Building $(BINARY_NAME)..."
	@mkdir -p bin
	@go build -o $(BINARY_PATH) .
	@echo "Built successfully: $(BINARY_PATH)"

# Clean build artifacts
clean:
	@echo "Cleaning..."
	@rm -f $(BINARY_PATH)
	@echo "Cleaned: $(BINARY_PATH)"

# Run tests
test:
	@echo "Running tests..."
	@go test ./...

# Install dependencies
deps:
	@echo "Installing dependencies..."
	@go mod download

# Run the binary
run: build
	@./$(BINARY_PATH)

# Default target
.DEFAULT_GOAL := build
