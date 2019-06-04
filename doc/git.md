# Git/GutHub Workflow

A quick guide to use git/github to contribute code to this MCCE project.

## Setup of git/github

### One time config - global settings:

Identity:

    git config --global user.name "Junjun Mao"
    git config --global user.email junjun.mao@gmail.com


Editor: 

    git config --global core.editor vi 

Save credential for 6 hours: 

    git config --global credential.helper 'cache --timeout=21600' 

Save credential permanently: 

    git config credential.helper store 
    git config --global credential.helper store 

Be aware, the credential is saved as clear text in file ```.git-credentials``` under project or in ```~/``` folder. 


## Common commands to work on your own repository/fork.
If you are developing code all by yourself, you don't have manage fork, branch and pull. This simplifies the work flow of using git and github.

This section also applies to project branch in the repository you already forked.

 
## Revisions - rewind the code history

## Solving conflict

## Collaboration in a project - Issue, Fork, Branch, Pull and Merge 

configure upstream
syncing a fork

push to the branch