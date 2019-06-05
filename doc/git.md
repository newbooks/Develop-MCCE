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

### To update and synchronize with github
Suppose you have a github repository "test", and cloned to office computer as "test". 

Now you added some files and updated some files. To check what are tracked:
```
git status
```
To add all files under current directory and subdirectories for tracking
```
git add .
```
To remove files/folders from tracking, there are two ways:
 * list the file or folder name in file ".gitignore" 
 * use ```git rm --cached filename``` to remove the file from tracking

Commit the changes:
```
git commit -a -m "remark of this commit"
```

Push the changes to github:
```
git push
```

Now you are back home and want to continue to work on home computer.

If not cloned yet, clone the project:
```
git clone https://github.com/newbooks/test.git
``` 
Otherwise, go to the project directory,
1 git pull
2 then edit code
3 git commit
4 git push 

Do 3 and 4 as often as you like.

Repeat 1 to 4 any time you switch a computer to continue to write code. This way, all computers are synchronized with github repository. 
 
## Revisions - rewind the code history

## Solving conflict

## Collaboration in a project - Issue, Fork, Branch, Pull and Merge 

configure upstream
syncing a fork

push to the branch
rebase to shorten the history