FROM ubuntu:22.04

# Install conda / mamba
RUN apt-get update && apt-get install -y wget

RUN CONDA="Mambaforge-pypy3-Linux-x86_64.sh" && \
    wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/$CONDA && \
    chmod +x $CONDA && \
    ./$CONDA -b -p /mambaforge-pypy3 && \
    rm -f $CONDA
ENV PATH /mambaforge-pypy3/bin:$PATH

# Install wic
RUN apt-get install -y git

COPY . /workflow-inference-compiler
WORKDIR /workflow-inference-compiler

#RUN conda create --name wic
#RUN conda activate wic
# The above command prints
# CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
# It still prints that even if we run `conda init bash` first.
# But this is a Docker image; we don't necessarily need to additionally isolate
# wic within a conda environment. Let's just install it globally!
RUN mamba env update --name base --file conda_devtools.yml
RUN pip install -e ".[all]"

ADD Dockerfile .
