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
If you are developing code all by yourself, you don't have to manage fork, branch and pull. The following work flow 
is for a single-developer project, which can be your original project or a forked project.

### To start a git project that synchronizes with github
* Create a repository (or fork a repository) on your github account <br>
  For example: Create a repository "test" under my account. <br>
  Login github -> https://github.com/new  -> give it a name "test"
* Get the remote url of this repository <br>
  Clone or download -> copy the url (https://github.com/newbooks/test.git in my case)
* Clone the project to local computer <br>
  got the directory that will contain the project, and clone the project <br>
  ```git clone https://github.com/newbooks/test.git```

*If you get an error message like this "WARNING: gnome-keyring:: couldn't connect to: ...", 
do the following under bash terminal:*
```    
unset GNOME_KEYRING_CONTROL
```

### To see what are tracked, add and remove files tracked

    git status
    git ls-files
    git add .
    git add file
    git rm --cached file

You can also ignore certain files and directories by putting file names and directories names in file .gitignore

### Working in a git directory
Commit the changes:
```
git commit -a -m "remark of this commit"
```

Push the changes to github:
```
git push
```

### Working on multiple computers to program the same project

If another computer doesn't have local repository directory yet, clone the project:
```
git clone https://github.com/newbooks/test.git
``` 
Otherwise, go to the project directory,
 1. git pull
 2. then edit code
 3. git commit
 4. git push 

Each commit saves a snapshot of the code and can be reverted in future.

With pull at the beginning and and push at the end, the changes are updated on github and all computers are 
synchronized. 
 
## Revisions - rewind the code history
### View history
One line short report: 

    git log --pretty=oneline 
    git log --pretty=format:"%h - %an, %ar : %s" 

Recent only: 

    git log --since=2weeks 
    git log --since="2015-1-15" 
    git log -10 

### Undo things
Discard changes on a file in working directory: 

    git checkout -- filename 
    git checkout <commit hash> path/to/file 

Revert to an old version. This command actually creates a new commit to image a previous state while keeping all the history. 
    git revert <commit hash>

You can use the first several characters of the commit hash, as long as it is enough to identify the commit.

### Remote repository
Show remote: 

    git remote 
    git remote -v 
    git remote show 
    git remote show <remote name>
    
Add a local directory to remote repository:
    git remote add origin https://github.com/newbooks/test.git
This command takes two arguments. The first is the remote name and the second is github repository URL. The remote 
name is a short label of remote repository. One can have a local git directory to be associated to more than one 
remote repository.

### Tagging
Tags are to mark release points. They are optional. To list tags: 
    git tag 
    git tag -l 'v8.5*' 
    git tag show 

Create light weight tag: 
    git tag v1.4

Create annotated tag with -a and optional -m 
    git tag -a v1.4 -m 'my version 1.4' 

Tag a commit: 

    git tag -a v1.2 9fceb02 

Push tag to remote server 

    git push origin [tagname] 


## Collaboration in a project - Issue, Fork, Branch, Pull and Merge 
When writing a program together, we usually work on different parts, or different aspects or the same program. One goal is to protect the master repository while allowing maximum flexibility for developers to experiment their implementation. So instead of editing the code in master repository, we work in forked repositories.

```
Upstream (master copy) -> fork -> issue -> branch -> coding -> pull request -> approval 
         |                                  |__________|                           |
         |________________________close issue, delete branch_______________________| 
                                            
```                    

Branch commands:

    git branch                          (show branches)
    git checkout -b <branch name>       (create and switch to the branch)
    git branch <branch name>            (switch branch)
    git branch -d <branch name>         (delete a branch)

configure upstream
syncing a fork

push to the branch
rebase to shorten the history