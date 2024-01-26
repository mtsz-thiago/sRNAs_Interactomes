FROM "mcr.microsoft.com/devcontainers/base:bullseye"

# Install NVIDIA container toolkit
RUN curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

RUN apt update && apt install -y \
    curl \
    git-flow \
    # git \
    python3 \
    python3-pip \
    python3-venv \
    graphviz \
    openjdk-11-jdk \
    nvidia-container-toolkit \
    && apt clean -y && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io > nextflow

RUN mv nextflow /usr/local/bin/
RUN chmod +x /usr/local/bin/nextflow
RUN nextflow
