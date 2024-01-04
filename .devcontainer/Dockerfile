FROM "mcr.microsoft.com/devcontainers/base:bullseye"

RUN apt update && apt install -y \
    curl \
    git-flow \
    # git \
    # python3 \
    # python3-pip \
    # python3-venv \
    openjdk-11-jdk \
    && apt clean -y && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io > nextflow

RUN mv nextflow /usr/local/bin/
RUN chmod +x /usr/local/bin/nextflow