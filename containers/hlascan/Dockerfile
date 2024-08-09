FROM ubuntu:24.10

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    wget \
    zip \
    build-essential \
    gcc \
    make \
    procps \
    samtools && \
    apt-get clean

RUN wget -O /usr/bin/hla_scan_r_v2.1.4 https://github.com/SyntekabioTools/HLAscan/releases/download/v2.1.4/hla_scan_r_v2.1.4
RUN mkdir -p /tmp
RUN wget -O /tmp/dataset.zip https://github.com/SyntekabioTools/HLAscan/releases/download/v2.0.0/dataset.zip
RUN unzip /tmp/dataset.zip
RUN chmod +x /usr/bin/hla_scan_r_v2.1.4

