FROM biocontainers/samtools:v1.9-4-deb_cv1

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    wget \
    zip \
    procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN wget -O /usr/bin/hla_scan_r_v2.1.4 https://github.com/SyntekabioTools/HLAscan/releases/download/v2.1.4/hla_scan_r_v2.1.4
RUN mkdir -p /tmp
RUN wget -O /tmp/dataset.zip https://github.com/SyntekabioTools/HLAscan/releases/download/v2.0.0/dataset.zip
RUN unzip /tmp/dataset.zip

