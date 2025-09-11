# Multi-stage build for Chromap with overflow system
# Stage 1: Build environment
FROM debian:trixie-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    make \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy source code
COPY . .

# Build Chromap with overflow system enabled
RUN make clean && make NEW_OVERFLOW=1 -j$(nproc)

# Verify the binary was built successfully
RUN ls -la chromap && ./chromap --help | head -5

# Stage 2: Runtime environment
FROM debian:trixie-slim AS runtime

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    zlib1g \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Create dedicated temp directory with proper permissions
RUN mkdir -p /chromap-temp && chmod 1777 /chromap-temp

# Set Chromap-specific temp directory environment variable
ENV CHROMAP_TEMP_DIR=/chromap-temp

# Copy the compiled binary from builder stage
COPY --from=builder /build/chromap /usr/local/bin/chromap

# Make binary executable
RUN chmod +x /usr/local/bin/chromap

# Create a non-root user for security
RUN useradd -r -s /bin/false chromap

# Set working directory
WORKDIR /data

# Create volume mount point for temp directory (optional)
VOLUME ["/chromap-temp"]

# Test that the binary works
RUN chromap --help | head -5

# Default command
CMD ["chromap", "--help"]

# Stage 3: Binary export stage (for extracting the binary)
FROM scratch AS export
COPY --from=builder /build/chromap /chromap

# Default stage is runtime (this must be last for docker build without --target)
FROM runtime
