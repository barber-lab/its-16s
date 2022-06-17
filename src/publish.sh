#!/usr/bin/env bash

#
# Publish results of an analysis
#

echo Saving results:
find results
rm results.zip
zip results.zip results

# get tag
release_tag=$CI_COMMIT_TAG
test -z $release_tag && \
	release_tag=$CI_COMMIT_SHORT_SHA && \
	git tag -a $release_tag -m "Autotag $(date)"

# Push results
echo push results
ls -l results
git add -A results
git	-c "user.name=Bot" \
	-c "user.email=bot@example.com" \
	commit -m "Update results: $(date)" 
git \
	-c http.sslverify="false" push -u --follow-tags \
	https://gitlab-ci-token:${GITLAB_ACCESS_TOKEN}@${CI_SERVER_HOST}/${CI_PROJECT_PATH}.git HEAD:master || true

# create result
echo create release
git_commit_with_results=$(git rev-parse HEAD)

cat << EOF > release.json
	{
		"name": "Results of analysis job #$CI_JOB_ID",
		"tag_name": "$release_tag",
		"description": "- CI/CD: Job [#$CI_JOB_ID]($CI_JOB_URL) in pipeline [#$CI_PIPELINE_ID]($CI_PIPELINE_URL)\n - Release date: $(date)\n - Docker image: [$CI_REGISTRY_IMAGE:$release_tag](https://sbi-git.hki-jena.de/$CI_PROJECT_PATH/container_registry)\n - Source code at analysis start: [Git commit $CI_COMMIT_SHA]($CI_SERVER_URL/$CI_PROJECT_PATH/-/tree/$CI_COMMIT_SHA)\n - Source code at analysis end with results: [Git commit $git_commit_with_results]($CI_SERVER_URL/$CI_PROJECT_PATH/-/tree/$git_commit_with_results)",
		"assets": {
			"links": [
				{
					"name": "results.zip",
					"url": "$CI_SERVER_URL/$CI_PROJECT_PATH/-/jobs/$CI_JOB_ID/artifacts/raw/results.zip",
					"link_type":"other"
				},
				{
					"name": "report.html",
					"url": "$CI_SERVER_URL/$CI_PROJECT_PATH/-/jobs/$CI_JOB_ID/artifacts/raw/results/report.html",
					"link_type":"other"
				},
				{
					"name": "manuscript.html",
					"url": "$CI_SERVER_URL/$CI_PROJECT_PATH/-/jobs/$CI_JOB_ID/artifacts/raw/results/manuscript.html",
					"link_type":"other"
				},
				{
					"name": "analyse.log",
					"url": "$CI_SERVER_URL/$CI_PROJECT_PATH/-/jobs/$CI_JOB_ID/artifacts/raw/log/analyse.log",
					"link_type":"other"
				}
			]
		}
	}
EOF
curl \
	-H 'Content-Type: application/json' \
	-H "PRIVATE-TOKEN: $GITLAB_ACCESS_TOKEN" \
	--data "@release.json" \
	-X POST \
	https://sbi-git.hki-jena.de/api/v4/projects/$CI_PROJECT_ID/releases