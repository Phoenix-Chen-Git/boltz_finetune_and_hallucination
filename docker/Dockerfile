# Use NVIDIA CUDA base image
FROM nvidia/cuda:12.1.1-devel-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV BOLTZ_CACHE=/root/.boltz

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-dev \
    python3-pip \
    git \
    wget \
    curl \
    aria2 \
    && rm -rf /var/lib/apt/lists/*

# Set python3.11 as default python
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --set python3 /usr/bin/python3.11 \
    && ln -s /usr/bin/python3 /usr/bin/python

# Upgrade pip and install build tools using python3.11 explicitly
RUN python3.11 -m pip install --upgrade pip setuptools wheel

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app

# Install boltz with CUDA support using the upgraded pip
RUN python3.11 -m pip install .[cuda]

# Install extra dependencies for processing scripts if needed
RUN if [ -f scripts/process/requirements.txt ]; then python3.11 -m pip install -r scripts/process/requirements.txt; fi

# Default command
ENTRYPOINT ["boltz"]
