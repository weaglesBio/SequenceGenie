FROM python:3.7

RUN apt-get -y update

# Setup paths:
ENV DIRPATH /sbc_ngs
WORKDIR $DIRPATH
COPY . .
ENV PYTHONPATH="$DIRPATH:$PYTHONPATH"

# Setup ARG variables:
ARG BWA_BIN="bwa-0.7.17.tar.bz2"
ARG BWA_VERSION="0.7.17"
ARG SAMTOOLS_BIN="samtools-1.6.tar.bz2"
ARG SAMTOOLS_VERSION="1.6"
ARG BCFTOOLS_BIN="bcftools-1.6.tar.bz2"
ARG BCFTOOLS_VERSION="1.6"
ARG MINIMAP2_BIN="minimap2-2.24.tar.bz2"
ARG MINIMAP2_VERSION="2.24"

# Install libraries:
RUN apt-get update && apt-get install -y --no-install-recommends \
	build-essential \
	ca-certificates \
	curl \
	libbz2-dev \
	liblzma-dev \
	libncurses5-dev \
	libncursesw5-dev \
	zlib1g-dev \
	&& rm -rf /var/lib/apt/lists/*

# Download and install minimap2:
RUN curl -fsSL https://github.com/lh3/minimap2/releases/download/v$MINIMAP2_VERSION/$MINIMAP2_BIN -o /opt/$MINIMAP2_BIN \
    && tar -jxvf /opt/$MINIMAP2_BIN -C /opt/ \
	&& cd /opt/minimap2-$MINIMAP2_VERSION \
	&& make \
	&& cp minimap2 /usr/local/bin \
	&& rm /opt/$MINIMAP2_BIN

#RUN curl -fsSL https://github.com/lh3/bwa/releases/download/v$BWA_VERSION/$BWA_BIN -o /opt/$BWA_BIN \
#	&& tar xjf /opt/$BWA_BIN -C /opt/ \
#	&& cd /opt/bwa-$BWA_VERSION \
#	&& pwd \
#	&& ls -l \
#	&& make \
#	&& mv bwa /usr/local/bin/ \
#	&& rm /opt/$BWA_BIN

# Download and install samtools:
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/$SAMTOOLS_BIN -o /opt/$SAMTOOLS_BIN \
	&& tar xjf /opt/$SAMTOOLS_BIN -C /opt/ \
	&& cd /opt/samtools-$SAMTOOLS_VERSION \
	&& make \
	&& make install \
	&& rm /opt/$SAMTOOLS_BIN

# Download and install bcftools:
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/$BCFTOOLS_BIN -o /opt/$BCFTOOLS_BIN \
	&& tar xjf /opt/$BCFTOOLS_BIN -C /opt/ \
	&& cd /opt/bcftools-$BCFTOOLS_VERSION \
	&& make \
	&& make install \
	&& rm /opt/$BCFTOOLS_BIN

# Install Python dependencies:
RUN pip install --upgrade pip \
  && pip install -r requirements.txt
  
# Set ENTRYPOINT:
ENTRYPOINT ["python", "-u", "app.py"]