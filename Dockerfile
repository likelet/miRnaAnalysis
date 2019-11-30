FROM nfcore/base
LABEL description="Docker image containing all requirements for ExomePipe pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sysucc-exomeseq-1.0dev/bin:$PATH