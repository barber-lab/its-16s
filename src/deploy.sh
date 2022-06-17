#!/usr/bin/env bash
echo Deployment started

export PATH=/usr/sbin:${PATH}
export USER=rstudio
groupadd docker
groupmod -g ${DOCKER_GROUP_ID} docker
usermod -u ${USER_ID} -g ${GROUP_ID} -d /analysis rstudio
cd /analysis
git fetch origin
git checkout master
git remote remove origin
# Set remote with access token having write access 
git remote add origin https://gitlab-ci-token:${GITLAB_ACCESS_TOKEN}@sbi-git.hki-jena.de/analyses/dloos/AX3-its-16s-stat.git
git fetch origin
git branch --set-upstream-to origin/master
# Prevent vars e.g. USER being forwared to RStudio to prevent log in confusion 
env | grep -E '^(CI|DOCKER|GIT)' | sed 's/=/="/' | sed 's/$/"/' >> /usr/local/lib/R/etc/Renviron
chown -R rstudio /analysis
chown rstudio /analysis/.git/config
chmod go-rwx /analysis/.git/config
/init
