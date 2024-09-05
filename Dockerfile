# Use CentOS 7 as the base image
FROM centos:7

# Import the CentOS 7 GPG key before installing packages
RUN rpm --import /etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7

# Fix CentOS mirror issues by using vault.centos.org
RUN sed -i 's|^mirrorlist=|#mirrorlist=|g' /etc/yum.repos.d/CentOS-Base.repo && \
    sed -i 's|^#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-Base.repo

# Install necessary tools and dependencies
RUN yum update -y && \
    yum install -y \
    epel-release \
    glibc-langpack-en \
    git \
    gcc \
    gcc-c++ \
    make \
    wget \
    zlib-devel \
    openssl-devel \
    curl-devel && \
    yum clean all

# Install a newer version of CMake
RUN wget https://cmake.org/files/v3.21/cmake-3.21.0-linux-x86_64.sh && \
    chmod +x cmake-3.21.0-linux-x86_64.sh && \
    ./cmake-3.21.0-linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm cmake-3.21.0-linux-x86_64.sh

# Set environment variables (use en_US.UTF-8 instead of C.UTF-8)
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# Clone UKFTractography from GitHub
WORKDIR /opt
RUN git clone https://github.com/pnlbwh/ukftractography.git

# Build UKFTractography
WORKDIR /opt/ukftractography
RUN mkdir build && cd build && \
    cmake .. && \
    make && \
    make test

# Clean up unnecessary files
RUN yum remove -y git gcc gcc-c++ make wget curl-devel && \
    yum clean all && \
    rm -rf /var/cache/yum

# Set the entrypoint for the container to run UKFTractography
# Adjust this path according to where the binary is actually built
ENTRYPOINT ["/opt/ukftractography/build/UKFTractography-build/UKFTractography/bin/UKFTractography"]

# Optional: Add default command or arguments if necessary
# CMD ["--help"]
