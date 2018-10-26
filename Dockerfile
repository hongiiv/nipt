FROM ubuntu:18.04

MAINTAINER cb.hong@ngenebio.com

LABEL www.ngenebio.com="NGeneBio" version="1.0.0" description="Docker NIPT"

USER root

RUN apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata && apt-get -qq -y install curl bzip2 && apt-get -qq -y install git r-base \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && conda config --append channels conda-forge \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

RUN conda install numpy scipy scikit-learn numexpr
RUN conda install -c bioconda bedtools

ENV PATH /opt/conda/bin:$PATH

ENV \
	INSTALL_PATH="/opt/ngenebio_ngs" \
	NGENEBIO_NGS_PATH="/opt/ngenebio_ngs/source" \
	MINICONDA_PATH="/opt/miniconda"
ENV \
	PATH="$NGENEBIO_NGS_PATH:$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH" \
	CONDA_DEFAULT_ENV=$MINICONDA_PATH \
	CONDA_PREFIX=$MINICONDA_PATH \
	JAVA_HOME=$MINICONDA_PATH

RUN mkdir -p $INSTALL_PATH
RUN mkdir -p $NGENEBIO_NGS_PATH
RUN mkdir -p $MINICONDA_PATH
RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

WORKDIR $INSTALL_PATH

COPY docker/install_ngenebio_ngs.sh $NGENEBIO_NGS_PATH/docker/
COPY requirements-minimal.txt $NGENEBIO_NGS_PATH/
RUN chmod +x $NGENEBIO_NGS_PATH/docker/install_ngenebio_ngs.sh 
COPY requirements-conda.txt requirements-conda-tests.txt requirements-py3.txt requirements.txt $NGENEBIO_NGS_PATH/
RUN $NGENEBIO_NGS_PATH/docker/install_ngenebio_ngs.sh minimal
RUN $NGENEBIO_NGS_PATH/docker/install_ngenebio_ngs.sh

COPY . $NGENEBIO_NGS_PATH

# Volume setup: make external tools and data available within container
VOLUME ["/ngenebio", "/user-data"]
ENV \
	NGENEBIO_NGS_DOCKER_DATA_PATH="/user-data" \
	GATK_PATH="/ngenebio"

RUN /bin/bash -c "set -e; check.py --version"

CMD ["/bin/bash"]
