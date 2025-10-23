# Stage 1: Build stage
FROM python:3.12-slim AS build

## Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
# setup python virtual environment
ENV PATH="/opt/venv/bin:$PATH"

# Set a working directory for building
WORKDIR /opt/proptimus

# Install system dependencies for building
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    curl \
    gfortran \
    git \
    libopenblas-dev \
    liblapack-dev \
    libeigen3-dev \
    libx11-dev \
    libglu1-mesa-dev \
    libxi-dev \
    libxrandr-dev \
    libxcursor-dev \
    libxtst-dev \
    zlib1g-dev

# Set up python virtual environment
RUN python3 -m venv /opt/venv

# Install xtb 6.6.1
RUN curl -L https://github.com/grimme-lab/xtb/releases/download/v6.6.1/xtb-6.6.1-source.tar.xz | tar xJ \
    && cd xtb-6.6.1 \
    && mkdir build && cd build \
    && cmake .. \
    && make -j$(nproc) \
    && make install \
    && cd ../.. && rm -rf xtb-6.6.1

# Install Python dependencies
RUN pip install --no-cache-dir \
    biopython==1.85 \
    rdkit==2025.09.1 \
    tqdm==4.67.1

## Get sources
COPY raphan.py .
COPY examples examples

### Stage 2: Runtime stage
FROM python:3.12-slim AS runtime

## Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/proptimus:/home/user/.local/bin:${PATH}"
# Use prepared python environment
ENV PATH="/opt/venv/bin:$PATH"

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libopenblas-dev \
    gfortran \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## Copy artefacts from build container
# Python libs
COPY --from=build /opt/venv /opt/venv
# xtb
COPY --from=build /usr/local/bin/xtb /usr/local/bin/
COPY --from=build /usr/local/lib/libxtb.so* /usr/local/lib/
COPY --from=build /usr/local/share/xtb /usr/local/share/xtb
# PDBCharges
COPY --from=build /opt/proptimus /opt/proptimus

# Set working directory
WORKDIR /opt/proptimus

# Create a non-root user and change ownership
RUN useradd --create-home --shell /bin/bash user \
    && chown -R user:user /opt \
    && chmod u+x raphan.py

# Switch to the non-root user
USER user

# Set default command
CMD ["python3", "raphan.py"]
