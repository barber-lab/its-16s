# get working with banocc and ComplexHeatmap
FROM rocker/verse:4.0.2

ENV DEBIAN_FRONTEND noninteractive
ENV TZ Europe/Berlin

USER root
WORKDIR /build

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && \
	apt-get install -y \
		make \
		# R sf https://github.com/r-spatial/
		libudunits2-dev \
		libgdal-dev \
		libgeos-dev \
		libproj-dev \
		git \
		python3-pip \
		htop \
		curl \
		zip \
		unzip \
		docker.io \
		parallel

RUN pip3 install dvc pre-commit

USER root
WORKDIR /build

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.2-Linux-x86_64.sh && \
	bash Miniconda3-py38_4.8.2-Linux-x86_64.sh -b -p /miniconda3 && \
 	/miniconda3/bin/conda init && \
	rm Miniconda*.sh && \
	chmod a+rwx -R /miniconda3
ENV PATH="/miniconda3/bin:${PATH}"

# install conda packages
RUN conda init bash &&  \
	conda create -y -n fastspar -c bioconda -c conda-forge fastspar mkl

COPY src/install/ src/install
RUN Rscript src/install/install.R

USER rstudio
EXPOSE 8787
WORKDIR /analysis
CMD make analyse
