#!/usr/bin/make -f

# This is the main entry point for the analysis
# All results must be reproducible by calling the rule 'run'

SHELL := /bin/bash

USER_ID := $(shell id -u)
GROUP_ID := $(shell id -g)

CPUS := $(shell nproc)
CACHE_DIR := ${PWD}/../cache

GIT_PROJECT_NAMESPACE := $(shell echo ${CI_PROJECT_NAMESPACE})
GIT_PROJECT_NAME := $(shell echo ${CI_PROJECT_NAME})
GIT_REPOSITORY_URL := $(shell echo ${CI_REPOSITORY_URL})

DOCKER_NAME := ax3-its-16s-stat
DOCKER_CLEANUP_STAGE := deploy
# ongoing numbering of deploy instances
DOCKER_DEPLOY_NR := $(shell docker container ls -a -f name=^/$(DOCKER_NAME)_deploy$ | wc -l)
DOCKER_DEPLOY_PORT := $(shell comm -23 <(seq 49152 65535 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | shuf | head -n 1)
DOCKER_DEPLOY_NAME := $(DOCKER_NAME)-deploy-$(DOCKER_DEPLOY_NR)
DOCKER_ANALYSE_NAME := $(DOCKER_NAME)-analyse-$(DOCKER_DEPLOY_NR)
DOCKER_REGISTRY := $(shell echo ${CI_PROJECT_NAME})
DOCKER_REGISTRY_USER := $(shell echo ${DOCKER_REGISTRY_USER})
DOCKER_REGISTRY_PASSWORD := $(shell echo ${DOCKER_REGISTRY_PASSWORD})
DOCKER_NAMESPACE := $(shell echo ${GIT_PROJECT_NAMESPACE} | tr [:upper:] [:lower:])
DOCKER_REPOSITORY := $(shell echo $(DOCKER_REGISTRY)/$(DOCKER_NAMESPACE)/${GIT_PROJECT_NAME} | sed 's/\/\//\//g' | tr '[:upper:]' '[:lower:]')

DOCKER_TAG := latest
DOCKER_IMAGE := $(DOCKER_REPOSITORY):$(DOCKER_TAG)
DOCKER_GROUP_ID := $(shell getent group | grep -E ^docker | cut -d : -f 3)

define DOCKER_RUN_ARGS
	--volume ${PWD}:/analysis \
	--volume $(CACHE_DIR):/cache \
	--volume /sbidata:/sbidata \
	--volume /scratch:/scratch \
	--volume /var/run/docker.sock:/var/run/docker.sock \
	--hostname $(DOCKER_NAME) \
	--env-file <(env) \
	--memory=500G \
	--cpus=$(CPUS)
endef

# Run docker container. This will also run the analysis by default
.PHONY: run
run: build
	docker tag $(DOCKER_REPOSITORY):$(DOCKER_TAG) $(DOCKER_IMAGE) && \
	docker push $(DOCKER_IMAGE) && \
	docker run \
		--user $(USER_ID):$(GROUP_ID) \
		--name $(DOCKER_ANALYSE_NAME) \
		--volume ${PWD}/src/Rprofile/run.R:/usr/local/lib/R/etc/Rprofile.site \
		$(DOCKER_RUN_ARGS) $(DOCKER_IMAGE) \
		bash -c "cd /analysis && make analyse"

# Start docker container for interactive analysis
.PHONY: deploy
deploy: build
	docker run \
		--detach \
		--user root:root \
		-e PASSWORD=password \
		-p $(DOCKER_DEPLOY_PORT):8787 \
		-v ${PWD}/src/Rprofile/env.R:/usr/local/lib/R/etc/Rprofile.site \
		--name $(DOCKER_DEPLOY_NAME) \
		-e GITLAB_ACCESS_TOKEN=${GITLAB_ACCESS_TOKEN} \
		-e DOCKER_GROUP_ID=$(DOCKER_GROUP_ID) \
		-e USER_ID=$(USER_ID) \
		-e GROUP_ID=$(GROUP_ID) \
		$(DOCKER_RUN_ARGS) $(DOCKER_IMAGE) \
		sh /analysis/src/deploy.sh \
	echo SUCCESS: Interactive analysis container $(DOCKER_DEPLOY_NAME) started.
	echo Access: http://$(shell hostname):$(DOCKER_DEPLOY_PORT) username:rstudio password:password

.PHONY: build
build:
	docker login \
		-u $(DOCKER_REGISTRY_USER) -p $(DOCKER_REGISTRY_PASSWORD) \
		$(DOCKER_REGISTRY) || true && \
	docker pull $(DOCKER_REPOSITORY):$(DOCKER_TAG) || true && \
	docker build \
	 	--tag $(DOCKER_REPOSITORY):$(DOCKER_TAG) \
		--cache-from $(DOCKER_REPOSITORY):$(DOCKER_TAG) \
		--build-arg BUILDKIT_INLINE_CACHE=1 \
		${PWD} \
		2>&1 | tee log/docker_build.log && \
	docker push $(DOCKER_REPOSITORY):$(DOCKER_TAG)

# Main entry point for the analysis
.PHONY: analyse
analyse:
	@echo "Start anlysis at ${HOSTNAME}"
	ls results | xargs -i rm -rf results/{}
	rm -rf log/*
	set -o pipefail && \
	Rscript src/main.R 2>&1 | tee log/analyse.log && \
	bash src/publish.sh 2>&1 | tee log/publish.log


.PHONY: cleanup
cleanup:
	docker rm -f $(shell docker ps -a -q --filter="name=$(DOCKER_CLEANUP_STAGE)-${CI_PIPELINE_ID}")
	rm -rf ../$(DOCKER_CLEANUP_STAGE)-${CI_PIPELINE_ID}*
